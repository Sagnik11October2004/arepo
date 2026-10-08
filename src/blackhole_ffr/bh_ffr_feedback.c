#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Iteration 7: bursty, exactly bipolar wind feedback.
 *
 * The full gas mass inside R_fb sets the binding-energy threshold, while only
 * hydro-active gas cells are eligible injection targets. Wind targets are
 * selected from the configured bicone, with a resolution fallback to larger
 * hemispheric receiver spheres. Each selected gas cell is kicked locally
 * outward along its BH-centered radial direction using mass-weighted top-hat
 * weights. The returned unresolved wind mass is added comoving with the BH
 * before the mechanical impulse is applied.
 *
 * Finite-resolution sampling can leave a small net gas impulse, so the BH
 * receives the equal-and-opposite recoil. The impulse amplitude is solved
 * from the exact combined gas+BH kinetic-energy quadratic, conserving total
 * momentum and injecting exactly E_packet up to MPI roundoff.
 *
 * The wind channel remains independent of the narrow jet channel. Iteration
 * 8 applies jet packets afterward using the gas state left by any wind events.
 */

#define BH_FFR_WIND_FALLBACK_RADIUS_COUNT 3

static const double BHFFRWindFallbackRadiusFactors[BH_FFR_WIND_FALLBACK_RADIUS_COUNT] =
    {64.0, 96.0, 128.0};

enum bh_ffr_feedback_pass
{
  BH_FFR_FEEDBACK_STATS = 0,
  BH_FFR_FEEDBACK_MARK = 1,
  BH_FFR_FEEDBACK_CONFLICT = 2,
  BH_FFR_FEEDBACK_INJECT = 3,
  BH_FFR_FEEDBACK_WAKE = 4
};

struct bh_ffr_feedback_event
{
  int Candidate;
  int Fire;
  MyDouble ThresholdEnergy;
  MyDouble PacketEnergy;
  MyDouble PacketMass;
  MyDouble PacketMomentum;
  MyDouble Q;
  MyDouble LobeMass[2];
  MyDouble LobeUtherm[2];
  long long LobeCount[2];
  MyDouble EnergyBefore;
  MyDouble MassBefore;
  MyDouble MomentumBefore;
  MyDouble RadiusProper;
  MyFloat RadiusCoordinate;
  MyFloat CosCone;
  int UsedHemisphereFallback;
  int RadiusIndex;
  int Underresolved;
};

static struct bh_ffr_feedback_event *FeedbackEvents;
static int *FeedbackPacketCount;
static MyIDType *FeedbackWinnerID;
static int FeedbackNTargets;
static int FeedbackPass;
static int FeedbackConflictRound;

typedef struct
{
  MyDouble Pos[3];
  MyFloat Axis[3];
  MyFloat BHVel[3];
  MyFloat Radius;
  MyDouble AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT];
  MyDouble FallbackRadius[BH_FFR_WIND_FALLBACK_RADIUS_COUNT];
  MyFloat CosCone;
  MyIDType BHID;
  int Candidate;
  int Fire;
  MyDouble PacketMass;
  MyDouble Q;
  MyDouble LobeMass[2];
  MyDouble LobeUtherm[2];
  int WakeTimeBin;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  /* Transient adaptive MACER wind-radius statistics. */
  MyDouble AdaptiveEnclosedMass[BH_FFR_ADAPTIVE_RADIUS_COUNT];
  MyDouble AdaptiveLobeMass[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  MyDouble AdaptiveLobeRadialMomentum[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  MyDouble AdaptiveLobeRadialDirectionMass[BH_FFR_ADAPTIVE_RADIUS_COUNT][2][3];
  MyDouble AdaptiveLobeThermalMassWeighted[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  MyDouble AdaptiveLobeTotalMass[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  long long AdaptiveLobeTotalCount[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  long long AdaptiveLobeCount[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];

  MyDouble FallbackEnclosedMass[BH_FFR_WIND_FALLBACK_RADIUS_COUNT];
  MyDouble FallbackHemiMass[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];
  MyDouble FallbackHemiRadialMomentum[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];
  MyDouble FallbackHemiRadialDirectionMass[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2][3];
  MyDouble FallbackHemiThermalMassWeighted[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];
  MyDouble FallbackHemiTotalMass[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];
  long long FallbackHemiTotalCount[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];
  long long FallbackHemiCount[BH_FFR_WIND_FALLBACK_RADIUS_COUNT][2];

  int Conflict;

  MyDouble ReturnedMass;
  MyDouble KickEnergy;
  MyDouble KickMomentum[3];
} data_out;

static data_out *FeedbackResults;
static data_out *DataResult, *DataOut;

static int bh_ffr_feedback_gas_is_active(int j)
{
  return bh_ffr_feedback_gas_is_hydro_active_now(j);
}

static int bh_ffr_feedback_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= FeedbackNTargets)
    terminate("BH_FFR: feedback target=%d outside [0,%d) in %s", target, FeedbackNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active feedback particle=%d in %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact feedback state for particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: feedback particle ID=%llu is not drifted to Ti_Current=%lld in %s", (unsigned long long)P[p].ID,
              (long long)All.Ti_Current, where);

  return p;
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_feedback_particle_from_target(target, "particle2in");
  const int b = P[p].BHDataIndex;
  const struct bh_ffr_feedback_event *ev = &FeedbackEvents[target];

  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[p].Pos[k];
      in->Axis[k] = BHP[b].DiscDir[k];
      in->BHVel[k] = P[p].Vel[k];
    }

  in->Radius = 0.0;
  for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
    in->AdaptiveRadius[k] = 0.0;
  for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
    in->FallbackRadius[k] = 0.0;

  if(FeedbackPass == BH_FFR_FEEDBACK_STATS)
    {
      struct bh_ffr_discrete_radius_grid grid, fallback_grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      bh_ffr_build_discrete_radius_grid(bh_ffr_effective_softening_proper(p),
                                        BHFFRWindFallbackRadiusFactors,
                                        BH_FFR_WIND_FALLBACK_RADIUS_COUNT,
                                        &fallback_grid);
      if(grid.Count != BH_FFR_ADAPTIVE_RADIUS_COUNT)
        terminate("BH_FFR: unexpected wind adaptive-radius count=%d", grid.Count);
      if(fallback_grid.Count != BH_FFR_WIND_FALLBACK_RADIUS_COUNT)
        terminate("BH_FFR: unexpected wind fallback-radius count=%d", fallback_grid.Count);

      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0))
        terminate("BH_FFR: invalid scale factor=%g in wind adaptive-radius setup", a);

      for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
        in->AdaptiveRadius[k] = grid.Radius[k] / a;
      for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
        in->FallbackRadius[k] = fallback_grid.Radius[k] / a;

      in->Radius = in->FallbackRadius[BH_FFR_WIND_FALLBACK_RADIUS_COUNT - 1];
    }
  else if(FeedbackPass == BH_FFR_FEEDBACK_WAKE)
    {
      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0))
        terminate("BH_FFR: invalid scale factor=%g in feedback sync guard", a);

      const double eps = bh_ffr_effective_softening_proper(p);
      const double rguard =
          BHFFRWindFallbackRadiusFactors[BH_FFR_WIND_FALLBACK_RADIUS_COUNT - 1] * eps;
      if(!isfinite(rguard) || !(rguard > 0))
        terminate("BH_FFR: invalid 128-epsilon feedback sync radius=%g ID=%llu",
                  rguard, (unsigned long long)P[p].ID);

      in->Radius = rguard / a;
    }
  else
    {
      if(!isfinite(ev->RadiusCoordinate) || !(ev->RadiusCoordinate > 0))
        terminate("BH_FFR: invalid selected wind coordinate radius=%g for ID=%llu",
                  (double)ev->RadiusCoordinate, (unsigned long long)P[p].ID);
      in->Radius = ev->RadiusCoordinate;
    }

  /* Statistics use the physical wind cone; injection uses the selected cone
   * or hemisphere; the shared synchronization guard is always a full sphere. */
  in->CosCone =
      FeedbackPass == BH_FFR_FEEDBACK_STATS
          ? cos(All.BHWindConeAngleDeg * M_PI / 180.0)
          : (FeedbackPass == BH_FFR_FEEDBACK_WAKE ? 0.0 : ev->CosCone);
  in->BHID = P[p].ID;
  in->WakeTimeBin =
      FeedbackPass == BH_FFR_FEEDBACK_WAKE
          ? bh_ffr_feedback_sync_target_hydro_timebin(p)
          : -1;
  in->Candidate = ev->Candidate;
  in->Fire = ev->Fire;
  in->PacketMass = ev->PacketMass;
  in->Q = ev->Q;
  in->LobeMass[0] = ev->LobeMass[0];
  in->LobeMass[1] = ev->LobeMass[1];
  in->LobeUtherm[0] = ev->LobeUtherm[0];
  in->LobeUtherm[1] = ev->LobeUtherm[1];
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= FeedbackNTargets)
    terminate("BH_FFR: feedback result target=%d outside [0,%d)", target, FeedbackNTargets);

  data_out *res = &FeedbackResults[target];

  if(mode == MODE_LOCAL_PARTICLES)
    *res = *out;
  else
    {
      for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
        {
          res->FallbackEnclosedMass[k] += out->FallbackEnclosedMass[k];
          for(int l = 0; l < 2; l++)
            {
              res->FallbackHemiMass[k][l] += out->FallbackHemiMass[k][l];
              res->FallbackHemiRadialMomentum[k][l] += out->FallbackHemiRadialMomentum[k][l];
              for(int q = 0; q < 3; q++)
                res->FallbackHemiRadialDirectionMass[k][l][q] +=
                    out->FallbackHemiRadialDirectionMass[k][l][q];
              res->FallbackHemiThermalMassWeighted[k][l] += out->FallbackHemiThermalMassWeighted[k][l];
              res->FallbackHemiTotalMass[k][l] += out->FallbackHemiTotalMass[k][l];
              res->FallbackHemiTotalCount[k][l] += out->FallbackHemiTotalCount[k][l];
              res->FallbackHemiCount[k][l] += out->FallbackHemiCount[k][l];
            }
        }
      for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
        {
          res->AdaptiveEnclosedMass[k] += out->AdaptiveEnclosedMass[k];
          for(int l = 0; l < 2; l++)
            {
              res->AdaptiveLobeMass[k][l] += out->AdaptiveLobeMass[k][l];
              res->AdaptiveLobeRadialMomentum[k][l] +=
                  out->AdaptiveLobeRadialMomentum[k][l];
              for(int q = 0; q < 3; q++)
                res->AdaptiveLobeRadialDirectionMass[k][l][q] +=
                    out->AdaptiveLobeRadialDirectionMass[k][l][q];
              res->AdaptiveLobeThermalMassWeighted[k][l] +=
                  out->AdaptiveLobeThermalMassWeighted[k][l];
              res->AdaptiveLobeTotalMass[k][l] +=
                  out->AdaptiveLobeTotalMass[k][l];
              res->AdaptiveLobeTotalCount[k][l] +=
                  out->AdaptiveLobeTotalCount[k][l];
              res->AdaptiveLobeCount[k][l] += out->AdaptiveLobeCount[k][l];
            }
        }
      if(out->Conflict)
        res->Conflict = 1;

      res->ReturnedMass += out->ReturnedMass;
      res->KickEnergy += out->KickEnergy;
      for(int k = 0; k < 3; k++)
        res->KickMomentum[k] += out->KickMomentum[k];
    }
}

#include "../utils/generic_comm_helpers2.h"

static int bh_ffr_feedback_evaluate(int target, int mode, int threadid);

static void kernel_local(void)
{
  const int threadid = get_thread_num();

  for(int task = 0; task < NTask; task++)
    Thread[threadid].Exportflag[task] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      const int target = NextParticle++;
      if(target >= FeedbackNTargets)
        break;

      bh_ffr_feedback_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_feedback_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static unsigned long long bh_ffr_feedback_priority(MyIDType id)
{
  /* Deterministic time-varying total order. This preserves conflict safety
   * without permanently starving the higher-ID BH in a persistent overlap. */
  unsigned long long x = (unsigned long long)id;
  x ^= (unsigned long long)All.Ti_Current + 0x9e3779b97f4a7c15ULL +
       (unsigned long long)(FeedbackConflictRound + 1) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static int bh_ffr_feedback_candidate_wins(MyIDType candidate, MyIDType incumbent)
{
  if(incumbent == 0)
    return 1;

  const unsigned long long pc = bh_ffr_feedback_priority(candidate);
  const unsigned long long pi = bh_ffr_feedback_priority(incumbent);
  return pc < pi || (pc == pi && candidate < incumbent);
}

static int bh_ffr_feedback_selected_lobe(const data_in *in, int j, double *r2_out)
{
  double xtmp, ytmp, ztmp;
  const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
  const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
  const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
  const double r2 = dx * dx + dy * dy + dz * dz;

  if(r2_out != NULL)
    *r2_out = r2;

  if(!(r2 > 0) || r2 > in->Radius * in->Radius)
    return -1;

  const double dot = dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
  const double cosine = fabs(dot) / sqrt(r2);

  if(cosine < in->CosCone)
    return -1;

  return dot >= 0 ? 0 : 1;
}

static int bh_ffr_feedback_evaluate(int target, int mode, int threadid)
{
  data_in local, *in;
  int numnodes, *firstnode;
  double xtmp, ytmp, ztmp;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  data_out out;
  memset(&out, 0, sizeof(out));

  if(FeedbackPass != BH_FFR_FEEDBACK_STATS)
    {
      if((FeedbackPass == BH_FFR_FEEDBACK_INJECT && !in->Fire) ||
         ((FeedbackPass == BH_FFR_FEEDBACK_MARK || FeedbackPass == BH_FFR_FEEDBACK_CONFLICT) && !in->Candidate))
        {
          if(mode == MODE_LOCAL_PARTICLES)
            out2particle(&out, target, MODE_LOCAL_PARTICLES);
          else
            DataResult[target] = out;
          return 0;
        }
    }

  const int nfound =
      ngb_treefind_variable_threads(in->Pos, in->Radius, target, mode, threadid, numnodes, firstnode);

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_FFR: invalid scale factor=%g in feedback evaluation", a);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_FFR: feedback gas index=%d outside NumGas=%d", j, NumGas);
      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      double r2;
      const int lobe = bh_ffr_feedback_selected_lobe(in, j, &r2);

      if(FeedbackPass == BH_FFR_FEEDBACK_WAKE)
        {
          if(lobe >= 0 && in->WakeTimeBin > 0)
            bh_ffr_feedback_request_wakeup(j, in->WakeTimeBin);
          continue;
        }

      if(FeedbackPass == BH_FFR_FEEDBACK_STATS)
        {
          const int active = bh_ffr_feedback_gas_is_active(j);
          int hemi_search = -1;

          if(r2 > 0 && r2 <= in->Radius * in->Radius)
            {
              const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
              const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
              const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
              const double dot = dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
              hemi_search = dot >= 0 ? 0 : 1;
            }

          /* One statistics pass builds every candidate sphere and bicone.
           * Thresholds use all enclosed gas; coupling quality is checked per
           * lobe using only hydro-active targets. */
          if(hemi_search >= 0 &&
             in->AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT - 1] > 0)
            {
              double radial_velocity = 0.0;
              double rhat[3] = {0.0, 0.0, 0.0};
              if(active)
                {
                  if(!isfinite(SphP[j].Utherm) || SphP[j].Utherm < 0)
                    terminate("BH_FFR: invalid target internal energy u=%g for gas ID=%llu",
                              (double)SphP[j].Utherm,
                              (unsigned long long)P[j].ID);
                  const double rinv = 1.0 / sqrt(r2);
                  rhat[0] = NEAREST_X(P[j].Pos[0] - in->Pos[0]) * rinv;
                  rhat[1] = NEAREST_Y(P[j].Pos[1] - in->Pos[1]) * rinv;
                  rhat[2] = NEAREST_Z(P[j].Pos[2] - in->Pos[2]) * rinv;
                  for(int q = 0; q < 3; q++)
                    radial_velocity += (P[j].Vel[q] / a) * rhat[q];
                }

              for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
                if(r2 <= in->AdaptiveRadius[k] * in->AdaptiveRadius[k])
                  {
                    out.AdaptiveEnclosedMass[k] += P[j].Mass;
                    if(lobe >= 0)
                      {
                        out.AdaptiveLobeTotalMass[k][lobe] += P[j].Mass;
                        out.AdaptiveLobeTotalCount[k][lobe]++;
                        if(active)
                          {
                            out.AdaptiveLobeCount[k][lobe]++;
                            out.AdaptiveLobeMass[k][lobe] += P[j].Mass;
                            out.AdaptiveLobeRadialMomentum[k][lobe] +=
                                P[j].Mass * radial_velocity;
                            for(int q = 0; q < 3; q++)
                              out.AdaptiveLobeRadialDirectionMass[k][lobe][q] +=
                                  P[j].Mass * rhat[q];
                            out.AdaptiveLobeThermalMassWeighted[k][lobe] +=
                                P[j].Mass * SphP[j].Utherm;
                          }
                      }
                  }

              for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
                if(in->FallbackRadius[k] > 0 && r2 <= in->FallbackRadius[k] * in->FallbackRadius[k])
                  {
                    out.FallbackEnclosedMass[k] += P[j].Mass;
                    out.FallbackHemiTotalMass[k][hemi_search] += P[j].Mass;
                    out.FallbackHemiTotalCount[k][hemi_search]++;
                    if(active)
                      {
                        out.FallbackHemiCount[k][hemi_search]++;
                        out.FallbackHemiMass[k][hemi_search] += P[j].Mass;
                        out.FallbackHemiRadialMomentum[k][hemi_search] += P[j].Mass * radial_velocity;
                        for(int q = 0; q < 3; q++)
                          out.FallbackHemiRadialDirectionMass[k][hemi_search][q] +=
                              P[j].Mass * rhat[q];
                        out.FallbackHemiThermalMassWeighted[k][hemi_search] += P[j].Mass * SphP[j].Utherm;
                      }
                  }
            }

          continue;
        }

      if(lobe < 0 || !bh_ffr_feedback_gas_is_active(j))
        continue;

      if(FeedbackPass == BH_FFR_FEEDBACK_MARK)
        {
          if(bh_ffr_feedback_candidate_wins(in->BHID, FeedbackWinnerID[j]))
            FeedbackWinnerID[j] = in->BHID;
          continue;
        }

      if(FeedbackPass == BH_FFR_FEEDBACK_CONFLICT)
        {
          if(FeedbackWinnerID[j] != in->BHID)
            out.Conflict = 1;
          continue;
        }

      if(FeedbackPass != BH_FFR_FEEDBACK_INJECT)
        terminate("BH_FFR: unknown feedback pass=%d", FeedbackPass);

      if(!(in->LobeMass[lobe] > 0))
        terminate("BH_FFR: non-positive injection lobe mass for BH ID=%llu", (unsigned long long)in->BHID);

      const double old_mass = P[j].Mass;
      const double weight = old_mass / in->LobeMass[lobe];
      const double dm = 0.5 * in->PacketMass * weight;

      if(!isfinite(weight) || !(weight > 0) || !isfinite(dm) || dm < 0)
        terminate("BH_FFR: invalid wind target weight=%g dm=%g for gas ID=%llu", weight, dm, (unsigned long long)P[j].ID);

      double bh_v2 = 0.0;
      for(int k = 0; k < 3; k++)
        bh_v2 += in->BHVel[k] * in->BHVel[k];

      P[j].Mass += dm;
      for(int k = 0; k < 3; k++)
        SphP[j].Momentum[k] += dm * in->BHVel[k];
      /* Return mass with the lobe-weighted ambient specific internal
       * energy; the mechanical packet remains a separate kinetic increment. */
      SphP[j].Energy += dm * (0.5 * bh_v2 + a * a * in->LobeUtherm[lobe]);

      const double rinv = 1.0 / sqrt(r2);
      const double rhat[3] = {
          NEAREST_X(P[j].Pos[0] - in->Pos[0]) * rinv,
          NEAREST_Y(P[j].Pos[1] - in->Pos[1]) * rinv,
          NEAREST_Z(P[j].Pos[2] - in->Pos[2]) * rinv};

      double dp_phys[3], v_phys[3];
      double dp2 = 0.0, vdotdp = 0.0;
      for(int k = 0; k < 3; k++)
        {
          v_phys[k] = (SphP[j].Momentum[k] / P[j].Mass) / a;
          dp_phys[k] = in->Q * weight * rhat[k];
          dp2 += dp_phys[k] * dp_phys[k];
          vdotdp += v_phys[k] * dp_phys[k];
        }

      const double dkin = vdotdp + 0.5 * dp2 / P[j].Mass;
      if(!isfinite(dkin))
        terminate("BH_FFR: non-finite wind kinetic increment for gas ID=%llu", (unsigned long long)P[j].ID);

      for(int k = 0; k < 3; k++)
        {
          SphP[j].Momentum[k] += a * dp_phys[k];
          out.KickMomentum[k] += dp_phys[k];
        }
      SphP[j].Energy += a * a * dkin;

      for(int k = 0; k < 3; k++)
        P[j].Vel[k] = SphP[j].Momentum[k] / P[j].Mass;

      out.ReturnedMass += dm;
      out.KickEnergy += dkin;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static void bh_ffr_feedback_comm_pass(int pass)
{
  FeedbackPass = pass;
  memset(FeedbackResults, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackResults));

  if(All.TotNumGas > 0)
    {
      generic_set_MaxNexport();
      generic_comm_pattern(FeedbackNTargets, kernel_local, kernel_imported);
    }
}

static double bh_ffr_feedback_packet_q(double A, double B, double energy)
{
  if(!isfinite(A) || !isfinite(B) || !(B > 0) || !isfinite(energy) || !(energy > 0))
    terminate("BH_FFR: invalid packet-root inputs A=%g B=%g E=%g", A, B, energy);

  const double disc = sqrt(A * A + 2.0 * B * energy);
  double q;

  if(A >= 0)
    q = 2.0 * energy / (disc + A);
  else
    q = (-A + disc) / B;

  if(!isfinite(q) || !(q > 0))
    terminate("BH_FFR: invalid positive packet root q=%g A=%g B=%g E=%g", q, A, B, energy);

  return q;
}

void bh_ffr_feedback_self_test(void)
{
  const double A = -0.37;
  const double B = 0.81;
  const double energy = 1.23;
  const double q = bh_ffr_feedback_packet_q(A, B, energy);
  const double got = A * q + 0.5 * B * q * q;

  if(fabs(got - energy) > 2.0e-13 * energy)
    terminate("BH_FFR: feedback self-test failed packet root got=%g expected=%g", got, energy);

  const double mplus = 2.0;
  const double mminus = 3.0;
  const double mpacket = 0.5;
  const double Bmass = 1.0 / (mplus + 0.5 * mpacket) + 1.0 / (mminus + 0.5 * mpacket);
  if(!isfinite(Bmass) || !(Bmass > 0))
    terminate("BH_FFR: feedback self-test failed mass-weighted B coefficient");

  /* Regression for PDF eq. (77): returned wind mass advects the lobe-weighted
   * ambient internal energy in addition to its BH-bulk kinetic energy. With no
   * mechanical kick and a common bulk velocity, this must preserve the
   * mass-weighted specific internal energy. The a^2 factors exercise AREPO's
   * cosmological conserved-energy convention directly. */
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_FFR: feedback self-test invalid scale factor=%g", a);

  const double m0 = 2.0, m1 = 3.0;
  const double u0 = 4.0, u1 = 10.0;
  const double m_before = m0 + m1;
  const double u_return = (m0 * u0 + m1 * u1) / m_before;
  const double v_bulk = 0.7;
  const double dm_return = 0.5;
  const double e_before = 0.5 * m_before * v_bulk * v_bulk + a * a * (m0 * u0 + m1 * u1);
  const double e_after =
      e_before + dm_return * (0.5 * v_bulk * v_bulk + a * a * u_return);
  const double m_after = m_before + dm_return;
  const double p_after = m_after * v_bulk;
  const double u_after =
      (e_after / m_after - 0.5 * (p_after / m_after) * (p_after / m_after)) / (a * a);

  if(!isfinite(u_after) || fabs(u_after - u_return) > 2.0e-13 * fmax(fabs(u_return), 1.0))
    terminate("BH_FFR: feedback self-test failed wind thermal return got=%g expected=%g", u_after, u_return);

  double previous_radius = 0.0;
  for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
    {
      const double r = BHFFRWindFallbackRadiusFactors[k];
      if(!isfinite(r) || !(r > 0) || (k > 0 && !(r > previous_radius)))
        terminate("BH_FFR: wind fallback-radius self-test failed k=%d R/eps=%g", k, r);
      previous_radius = r;
    }

  if(BHFFRWindFallbackRadiusFactors[0] != 64.0 ||
     BHFFRWindFallbackRadiusFactors[BH_FFR_WIND_FALLBACK_RADIUS_COUNT - 1] != 128.0)
    terminate("BH_FFR: wind fallback-radius endpoints are not 64 and 128 epsilon");

  const double Arad = -0.21;
  const double Bgas = 0.73;
  const double s2 = 0.14;
  const double mbh = 11.0;
  const double Erad = 0.91;
  const double Brad = Bgas + s2 / mbh;
  const double qrad = bh_ffr_feedback_packet_q(Arad, Brad, Erad);
  const double Erad_got = Arad * qrad + 0.5 * Brad * qrad * qrad;
  if(fabs(Erad_got - Erad) > 2.0e-13 * Erad)
    terminate("BH_FFR: radial-wind recoil self-test failed got=%g expected=%g",
              Erad_got, Erad);
}

static int bh_ffr_feedback_radius_geometry_ok(const data_out *res, int k)
{
  return res->AdaptiveLobeTotalCount[k][0] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->AdaptiveLobeTotalCount[k][1] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->AdaptiveLobeTotalMass[k][0] > 0 &&
         res->AdaptiveLobeTotalMass[k][1] > 0;
}

static int bh_ffr_feedback_radius_active_ok(const data_out *res, int k,
                                            double *fplus, double *fminus)
{
  const double mtot_plus = res->AdaptiveLobeTotalMass[k][0];
  const double mtot_minus = res->AdaptiveLobeTotalMass[k][1];
  *fplus = mtot_plus > 0 ? res->AdaptiveLobeMass[k][0] / mtot_plus : 0.0;
  *fminus = mtot_minus > 0 ? res->AdaptiveLobeMass[k][1] / mtot_minus : 0.0;

  return res->AdaptiveLobeCount[k][0] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->AdaptiveLobeCount[k][1] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         *fplus >= All.BHMinActiveTargetMassFrac &&
         *fminus >= All.BHMinActiveTargetMassFrac;
}

static int bh_ffr_feedback_fallback_geometry_ok(const data_out *res, int k)
{
  return res->FallbackHemiTotalCount[k][0] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->FallbackHemiTotalCount[k][1] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->FallbackHemiTotalMass[k][0] > 0 &&
         res->FallbackHemiTotalMass[k][1] > 0;
}

static int bh_ffr_feedback_fallback_active_ok(const data_out *res, int k,
                                               double *fplus, double *fminus)
{
  const double mtot_plus = res->FallbackHemiTotalMass[k][0];
  const double mtot_minus = res->FallbackHemiTotalMass[k][1];
  *fplus = mtot_plus > 0 ? res->FallbackHemiMass[k][0] / mtot_plus : 0.0;
  *fminus = mtot_minus > 0 ? res->FallbackHemiMass[k][1] / mtot_minus : 0.0;

  return res->FallbackHemiCount[k][0] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->FallbackHemiCount[k][1] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         *fplus >= All.BHMinActiveTargetMassFrac &&
         *fminus >= All.BHMinActiveTargetMassFrac;
}

static int bh_ffr_feedback_select_live_geometry(const data_out *res,
                                                int *fallback,
                                                int *underresolved,
                                                double *fplus,
                                                double *fminus)
{
  *fallback = 0;
  *underresolved = 0;
  *fplus = 0.0;
  *fminus = 0.0;

  for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
    {
      if(bh_ffr_feedback_radius_geometry_ok(res, k))
        {
          bh_ffr_feedback_radius_active_ok(res, k, fplus, fminus);
          return k;
        }
    }

  *fallback = 1;
  for(int k = 0; k < BH_FFR_WIND_FALLBACK_RADIUS_COUNT; k++)
    {
      if(bh_ffr_feedback_fallback_geometry_ok(res, k))
        {
          bh_ffr_feedback_fallback_active_ok(res, k, fplus, fminus);
          return k;
        }
    }

  const int k = BH_FFR_WIND_FALLBACK_RADIUS_COUNT - 1;
  bh_ffr_feedback_fallback_active_ok(res, k, fplus, fminus);
  *underresolved = 1;
  return k;
}

static void bh_ffr_feedback_report_adaptive_radius_survey(void)
{
  if(All.BHBenchmarkFeedbackModel != BH_BENCHMARK_FEEDBACK_MACER || FeedbackConflictRound != 0)
    return;
  static unsigned long long survey_calls = 0;
  if((survey_calls++ % 64ULL) != 0ULL)
    return;

  for(int n = 0; n < FeedbackNTargets; n++)
    {
      const int p = bh_ffr_feedback_particle_from_target(n, "bh_ffr_feedback_report_adaptive_radius_survey");
      const data_out *res = &FeedbackResults[n];
      const double eps = bh_ffr_effective_softening_proper(p);
      struct bh_ffr_discrete_radius_grid grid, fallback_grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      bh_ffr_build_discrete_radius_grid(eps, BHFFRWindFallbackRadiusFactors,
                                        BH_FFR_WIND_FALLBACK_RADIUS_COUNT, &fallback_grid);

      for(int k = 0; k < grid.Count; k++)
        {
          double fp = 0.0, fm = 0.0;
          const int geometry_ok = bh_ffr_feedback_radius_geometry_ok(res, k);
          const int active_ok = bh_ffr_feedback_radius_active_ok(res, k, &fp, &fm);
          printf("BH_FFR: adaptive wind radius detail mode=live ID=%llu task=%d geometry=cone index=%d R=%g Reps=%g NtotPlus=%lld NtotMinus=%lld NactivePlus=%lld NactiveMinus=%lld fplus=%g fminus=%g geometryResolved=%d activeReady=%d\n",
                 (unsigned long long)P[p].ID, ThisTask, k, (double)grid.Radius[k],
                 (double)(grid.Radius[k] / eps),
                 res->AdaptiveLobeTotalCount[k][0], res->AdaptiveLobeTotalCount[k][1],
                 res->AdaptiveLobeCount[k][0], res->AdaptiveLobeCount[k][1],
                 fp, fm, geometry_ok, active_ok);
        }
      for(int k = 0; k < fallback_grid.Count; k++)
        {
          double fp = 0.0, fm = 0.0;
          const int geometry_ok = bh_ffr_feedback_fallback_geometry_ok(res, k);
          const int active_ok = bh_ffr_feedback_fallback_active_ok(res, k, &fp, &fm);
          printf("BH_FFR: adaptive wind radius detail mode=live ID=%llu task=%d geometry=hemisphere index=%d R=%g Reps=%g NtotPlus=%lld NtotMinus=%lld NactivePlus=%lld NactiveMinus=%lld fplus=%g fminus=%g geometryResolved=%d activeReady=%d\n",
                 (unsigned long long)P[p].ID, ThisTask, k, (double)fallback_grid.Radius[k],
                 (double)(fallback_grid.Radius[k] / eps),
                 res->FallbackHemiTotalCount[k][0], res->FallbackHemiTotalCount[k][1],
                 res->FallbackHemiCount[k][0], res->FallbackHemiCount[k][1],
                 fp, fm, geometry_ok, active_ok);
        }

      int fallback = 0, underresolved = 0;
      double fp = 0.0, fm = 0.0;
      const int selected = bh_ffr_feedback_select_live_geometry(res, &fallback, &underresolved, &fp, &fm);
      const double rsel = fallback ? fallback_grid.Radius[selected] : grid.Radius[selected];
      const long long ntp = fallback ? res->FallbackHemiTotalCount[selected][0] : res->AdaptiveLobeTotalCount[selected][0];
      const long long ntm = fallback ? res->FallbackHemiTotalCount[selected][1] : res->AdaptiveLobeTotalCount[selected][1];
      const long long nap = fallback ? res->FallbackHemiCount[selected][0] : res->AdaptiveLobeCount[selected][0];
      const long long nam = fallback ? res->FallbackHemiCount[selected][1] : res->AdaptiveLobeCount[selected][1];
      const int active_ready =
          !underresolved &&
          (fallback ? bh_ffr_feedback_fallback_active_ok(res, selected, &fp, &fm)
                    : bh_ffr_feedback_radius_active_ok(res, selected, &fp, &fm));
      printf("BH_FFR: adaptive wind radius mode=live ID=%llu task=%d eps=%g Rsearch=%g Rwind=%g index=%d NtotPlus=%lld NtotMinus=%lld NactivePlus=%lld NactiveMinus=%lld fplus=%g fminus=%g fallback=%d underresolved=%d activeReady=%d\n",
             (unsigned long long)P[p].ID, ThisTask, eps,
             (double)fallback_grid.Radius[fallback_grid.Count - 1], rsel, selected,
             ntp, ntm, nap, nam, fp, fm, fallback, underresolved, active_ready);
      fflush(stdout);
    }
}

static void bh_ffr_feedback_prepare_candidates(void)
{
  for(int n = 0; n < FeedbackNTargets; n++)
    {
      struct bh_ffr_feedback_event *ev = &FeedbackEvents[n];
      memset(ev, 0, sizeof(*ev));
      const int p = bh_ffr_feedback_particle_from_target(n, "bh_ffr_feedback_prepare_candidates");
      const int b = P[p].BHDataIndex;
      const data_out *res = &FeedbackResults[n];

      struct bh_ffr_discrete_radius_grid grid, fallback_grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      bh_ffr_build_discrete_radius_grid(bh_ffr_effective_softening_proper(p),
                                        BHFFRWindFallbackRadiusFactors,
                                        BH_FFR_WIND_FALLBACK_RADIUS_COUNT, &fallback_grid);

      int fallback = 0, underresolved = 0;
      double fplus = 0.0, fminus = 0.0;
      const int selected = bh_ffr_feedback_select_live_geometry(res, &fallback, &underresolved, &fplus, &fminus);
      ev->RadiusIndex = selected;
      ev->Underresolved = underresolved;
      ev->UsedHemisphereFallback = fallback;
      ev->RadiusProper = fallback ? fallback_grid.Radius[selected] : grid.Radius[selected];

      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0)) terminate("BH_FFR: invalid scale factor=%g in wind candidate setup", a);
      ev->RadiusCoordinate = ev->RadiusProper / a;
      ev->CosCone = fallback ? 0.0 : cos(All.BHWindConeAngleDeg * M_PI / 180.0);

      const double menc = fallback ? res->FallbackEnclosedMass[selected] : res->AdaptiveEnclosedMass[selected];
      const double sigma2 = BHP[b].SigmaDM * BHP[b].SigmaDM;
      double vbind2 = sigma2;
      if(All.BHUseCentralBindingTerm) vbind2 += 2.0 * All.G * bh_ffr_central_mass_code(p) / ev->RadiusProper;
      if(!isfinite(vbind2) || vbind2 < 0 || !isfinite(menc) || menc < 0)
        terminate("BH_FFR: invalid adaptive wind threshold environment Menc=%g R=%g vbind2=%g for ID=%llu", menc, ev->RadiusProper, vbind2, (unsigned long long)P[p].ID);

      BHP[b].WindThresholdEnergy = All.BHWindBurstFactor * 0.5 * menc * vbind2;

      const int active_ready =
          !underresolved &&
          (fallback ? bh_ffr_feedback_fallback_active_ok(res, selected, &fplus, &fminus)
                    : bh_ffr_feedback_radius_active_ok(res, selected, &fplus, &fminus));

      if(BHP[b].SigmaDM > 0)
        {
          printf("BH_FFR: binding threshold wind ID=%llu task=%d Rwind=%g index=%d fallback=%d Menc=%g sigmaDM=%g vbind2=%g central=%d EthWind=%g\n",
                 (unsigned long long)P[p].ID, ThisTask, ev->RadiusProper, selected, fallback,
                 menc, BHP[b].SigmaDM, vbind2, All.BHUseCentralBindingTerm, BHP[b].WindThresholdEnergy);
          fflush(stdout);
        }
      if(underresolved || !active_ready) continue;

      const double eth = BHP[b].WindThresholdEnergy;
      if(!(eth > 0) || !(BHP[b].WindEnergyBuffer >= eth)) continue;
      ev->EnergyBefore = BHP[b].WindEnergyBuffer;
      ev->MassBefore = BHP[b].WindMassBuffer;
      ev->MomentumBefore = BHP[b].WindMomentumBuffer;
      ev->ThresholdEnergy = eth;
      ev->PacketEnergy = eth;
      const double frac = ev->PacketEnergy / ev->EnergyBefore;
      ev->PacketMass = ev->MassBefore * frac;
      ev->PacketMomentum = ev->MomentumBefore * frac;

      for(int l = 0; l < 2; l++)
        {
          if(fallback)
            {
              ev->LobeMass[l] = res->FallbackHemiMass[selected][l];
              ev->LobeCount[l] = res->FallbackHemiCount[selected][l];
              ev->LobeUtherm[l] = res->FallbackHemiThermalMassWeighted[selected][l] / ev->LobeMass[l];
            }
          else
            {
              ev->LobeMass[l] = res->AdaptiveLobeMass[selected][l];
              ev->LobeCount[l] = res->AdaptiveLobeCount[selected][l];
              ev->LobeUtherm[l] = res->AdaptiveLobeThermalMassWeighted[selected][l] / ev->LobeMass[l];
            }
          if(!isfinite(ev->LobeUtherm[l]) || ev->LobeUtherm[l] < 0)
            terminate("BH_FFR: invalid adaptive wind lobe internal energy=%g for ID=%llu lobe=%d", ev->LobeUtherm[l], (unsigned long long)P[p].ID, l);
        }

      if(!isfinite(frac) || frac <= 0 || frac > 1.0 + 2.0e-12 || !isfinite(ev->PacketMass) || ev->PacketMass < 0 || ev->PacketMass > ev->MassBefore * (1.0 + 2.0e-12))
        terminate("BH_FFR: invalid wind packet fraction=%g dM=%g for ID=%llu", frac, ev->PacketMass, (unsigned long long)P[p].ID);

      const double half = 0.5 * ev->PacketMass;
      const double mplus = ev->LobeMass[0] + half, mminus = ev->LobeMass[1] + half;
      if(!isfinite(mplus) || !isfinite(mminus) || !(mplus > 0) || !(mminus > 0))
        terminate("BH_FFR: invalid selected wind lobe masses ID=%llu fallback=%d mplus=%g mminus=%g", (unsigned long long)P[p].ID, fallback, mplus, mminus);

      double vbh_phys[3], total_direction[3] = {0.0, 0.0, 0.0};
      for(int k = 0; k < 3; k++)
        vbh_phys[k] = P[p].Vel[k] / a;

      double A = 0.0;
      double B = 1.0 / mplus + 1.0 / mminus;

      for(int l = 0; l < 2; l++)
        {
          const double mlobe = ev->LobeMass[l];
          const double mpost = mlobe + half;
          const double radial_momentum =
              fallback ? res->FallbackHemiRadialMomentum[selected][l]
                       : res->AdaptiveLobeRadialMomentum[selected][l];

          double mean_rhat[3];
          double vbh_dot_rhat = 0.0;
          for(int k = 0; k < 3; k++)
            {
              const double direction_mass =
                  fallback ? res->FallbackHemiRadialDirectionMass[selected][l][k]
                           : res->AdaptiveLobeRadialDirectionMass[selected][l][k];
              mean_rhat[k] = direction_mass / mlobe;
              total_direction[k] += mean_rhat[k];
              vbh_dot_rhat += vbh_phys[k] * mean_rhat[k];
            }

          A += (radial_momentum + half * vbh_dot_rhat) / mpost;
        }

      const double bh_mass_after = P[p].Mass - ev->PacketMass;
      if(!isfinite(bh_mass_after) || !(bh_mass_after > 0))
        terminate("BH_FFR: invalid BH dynamical mass after radial wind packet ID=%llu Mafter=%g",
                  (unsigned long long)P[p].ID, bh_mass_after);

      double total_direction2 = 0.0, vbh_dot_total = 0.0;
      for(int k = 0; k < 3; k++)
        {
          total_direction2 += total_direction[k] * total_direction[k];
          vbh_dot_total += vbh_phys[k] * total_direction[k];
        }

      A -= vbh_dot_total;
      B += total_direction2 / bh_mass_after;

      ev->Q = bh_ffr_feedback_packet_q(A, B, ev->PacketEnergy);
      ev->Candidate = 1;
    }
}

static int bh_ffr_feedback_global_candidate_count(void)
{
  int local = 0;
  for(int n = 0; n < FeedbackNTargets; n++)
    if(FeedbackEvents[n].Candidate)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static int bh_ffr_feedback_global_fire_count(void)
{
  int local = 0;
  for(int n = 0; n < FeedbackNTargets; n++)
    if(FeedbackEvents[n].Fire)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static void bh_ffr_feedback_resolve_conflicts(void)
{
  const int candidates = bh_ffr_feedback_global_candidate_count();
  if(candidates <= 0)
    return;

  if(candidates == 1)
    {
      for(int n = 0; n < FeedbackNTargets; n++)
        FeedbackEvents[n].Fire = FeedbackEvents[n].Candidate;
      return;
    }

  memset(FeedbackWinnerID, 0, (NumGas > 0 ? NumGas : 1) * sizeof(*FeedbackWinnerID));
  bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_MARK);
  bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_CONFLICT);

  for(int n = 0; n < FeedbackNTargets; n++)
    FeedbackEvents[n].Fire = FeedbackEvents[n].Candidate && !FeedbackResults[n].Conflict;
}

static void bh_ffr_feedback_commit_packets(void)
{
  for(int n = 0; n < FeedbackNTargets; n++)
    {
      struct bh_ffr_feedback_event *ev = &FeedbackEvents[n];
      if(!ev->Fire)
        continue;

      const int p = bh_ffr_feedback_particle_from_target(n, "bh_ffr_feedback_commit_packets");
      const int b = P[p].BHDataIndex;
      const data_out *res = &FeedbackResults[n];

      const double mtol = 2.0e-8 * fmax(fabs(ev->PacketMass), 1.0e-30);
      if(fabs(res->ReturnedMass - ev->PacketMass) > mtol)
        terminate("BH_FFR: wind packet mass-return mismatch ID=%llu returned=%g requested=%g", (unsigned long long)P[p].ID,
                  res->ReturnedMass, ev->PacketMass);

      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0))
        terminate("BH_FFR: invalid scale factor=%g in radial wind commit", a);

      const double bh_mass_after =
          BHP[b].BHMass + BHP[b].ReservoirMass +
          (ev->MassBefore - ev->PacketMass);
      if(!isfinite(bh_mass_after) || !(bh_mass_after > 0))
        terminate("BH_FFR: invalid BH mass after radial wind return ID=%llu M=%g",
                  (unsigned long long)P[p].ID, bh_mass_after);

      double gas_p2 = 0.0, bh_recoil_energy = 0.0;
      for(int k = 0; k < 3; k++)
        {
          const double dp_bh = -res->KickMomentum[k];
          const double vbh_phys = P[p].Vel[k] / a;
          gas_p2 += res->KickMomentum[k] * res->KickMomentum[k];
          bh_recoil_energy +=
              vbh_phys * dp_bh +
              0.5 * dp_bh * dp_bh / bh_mass_after;
        }

      const double total_kick_energy = res->KickEnergy + bh_recoil_energy;
      const double etol = 3.0e-8 * fmax(fabs(ev->PacketEnergy), 1.0e-30);
      if(fabs(total_kick_energy - ev->PacketEnergy) > etol)
        terminate("BH_FFR: radial wind packet energy mismatch ID=%llu gas=%g recoil=%g total=%g requested=%g",
                  (unsigned long long)P[p].ID, res->KickEnergy,
                  bh_recoil_energy, total_kick_energy, ev->PacketEnergy);

      BHP[b].WindEnergyBuffer = ev->EnergyBefore - ev->PacketEnergy;
      BHP[b].WindMassBuffer = ev->MassBefore - ev->PacketMass;
      BHP[b].WindMomentumBuffer = ev->MomentumBefore - ev->PacketMomentum;

      if(BHP[b].WindEnergyBuffer < 0 && BHP[b].WindEnergyBuffer > -3.0e-12 * fmax(fabs(ev->EnergyBefore), 1.0e-30))
        BHP[b].WindEnergyBuffer = 0.0;
      if(BHP[b].WindMassBuffer < 0 && BHP[b].WindMassBuffer > -3.0e-12 * fmax(fabs(ev->MassBefore), 1.0e-30))
        BHP[b].WindMassBuffer = 0.0;
      if(BHP[b].WindMomentumBuffer < 0 && BHP[b].WindMomentumBuffer > -3.0e-12 * fmax(fabs(ev->MomentumBefore), 1.0e-30))
        BHP[b].WindMomentumBuffer = 0.0;

      if(BHP[b].WindEnergyBuffer < 0 || BHP[b].WindMassBuffer < 0 || BHP[b].WindMomentumBuffer < 0)
        terminate("BH_FFR: negative wind buffer after packet for ID=%llu", (unsigned long long)P[p].ID);

      P[p].Mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
      if(fabs(P[p].Mass - bh_mass_after) >
         3.0e-12 * fmax(fabs(bh_mass_after), 1.0e-30))
        terminate("BH_FFR: radial wind BH mass mismatch ID=%llu state=%g expected=%g",
                  (unsigned long long)P[p].ID, (double)P[p].Mass, bh_mass_after);

      for(int k = 0; k < 3; k++)
        P[p].Vel[k] += a * (-res->KickMomentum[k]) / P[p].Mass;

      FeedbackPacketCount[n]++;

      const double gas_pnorm = sqrt(gas_p2);
      const double recoil_fraction = gas_pnorm / fmax(fabs(ev->Q), 1.0e-30);
      const double eerg = ev->PacketEnergy * All.UnitEnergy_in_cgs / All.HubbleParam;
      const double mmsun = ev->PacketMass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
      printf("BH_FFR: wind packet fired ID=%llu task=%d E=%g erg dM=%g Msun q=%g Nplus=%lld Nminus=%lld dErel=%g pbal=0 recoilFrac=%g bhRecoilEfrac=%g\n",
             (unsigned long long)P[p].ID, ThisTask, eerg, mmsun, ev->Q,
             ev->LobeCount[0], ev->LobeCount[1],
             fabs(total_kick_energy - ev->PacketEnergy) / ev->PacketEnergy,
             recoil_fraction, bh_recoil_energy / ev->PacketEnergy);
      printf("BH_FFR: wind coupling ID=%llu task=%d kick=radial-outward Rwind=%g index=%d "
             "fallback=%d cosCone=%g Nplus=%lld Nminus=%lld\n",
             (unsigned long long)P[p].ID, ThisTask, ev->RadiusProper,
             ev->RadiusIndex, ev->UsedHemisphereFallback,
             (double)ev->CosCone, ev->LobeCount[0], ev->LobeCount[1]);
      fflush(stdout);
    }
}

void bh_ffr_inject_wind_feedback(void)
{
  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

#ifdef MHD
  terminate("BH_FFR: kinetic wind feedback is not enabled with MHD until the magnetic mass-return policy is validated");
#endif

  FeedbackNTargets = NumActiveBHFFR;
  FeedbackPacketCount =
      (int *)mymalloc("BHFFRFeedbackPacketCount", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(int));
  FeedbackEvents = (struct bh_ffr_feedback_event *)mymalloc(
      "BHFFRFeedbackEvents", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackEvents));
  FeedbackResults =
      (data_out *)mymalloc("BHFFRFeedbackResults", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackResults));
  FeedbackWinnerID =
      (MyIDType *)mymalloc("BHFFRFeedbackWinner", (NumGas > 0 ? NumGas : 1) * sizeof(*FeedbackWinnerID));

  memset(FeedbackPacketCount, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(int));
  memset(FeedbackEvents, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackEvents));

  int any_global_fire = 0;

  for(int round = 0; round < BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP; round++)
    {
      FeedbackConflictRound = round;
      bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_STATS);
      bh_ffr_feedback_report_adaptive_radius_survey();
      bh_ffr_feedback_prepare_candidates();

      /* The shared 128-epsilon synchronization guard is executed once
       * in the jet stage, after both the current wind and current jet
       * thresholds have been surveyed. */

      if(bh_ffr_feedback_global_candidate_count() <= 0)
        break;

      bh_ffr_feedback_resolve_conflicts();
      const int global_fire = bh_ffr_feedback_global_fire_count();
      if(global_fire <= 0)
        break;

      bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_INJECT);
      bh_ffr_feedback_commit_packets();
      any_global_fire += global_fire;

      /* Wind coupling uses lobe-averaged Utherm on the next round. Refresh
       * primitives now so repeated packets never use stale thermodynamics. */
      update_primitive_variables();
    }

  /* End-of-transaction drainage diagnostic. A residual ratio >= 1 means
   * at least one threshold quantum remains buffered. The per-step packet
   * ceiling deliberately leaves any excess for later hydro updates so a
   * recovered coupling geometry cannot dump a long backlog into one state. */
  for(int n = 0; n < FeedbackNTargets; n++)
    {
      const int p =
          bh_ffr_feedback_particle_from_target(
              n, "bh_ffr_inject_wind_feedback/drain");
      const int b = P[p].BHDataIndex;
      const double eth = BHP[b].WindThresholdEnergy;
      const double ratio =
          eth > 0 ? BHP[b].WindEnergyBuffer / eth : 0.0;
      const int backlog = eth > 0 && BHP[b].WindEnergyBuffer >= eth;
      const int cap_hit =
          FeedbackPacketCount[n] >= BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP;

      printf("BH_FFR: wind drain ID=%llu task=%d packets=%d maxPackets=%d "
             "buffer=%g Eth=%g backlogRatio=%g backlog=%d capHit=%d\n",
             (unsigned long long)P[p].ID, ThisTask,
             FeedbackPacketCount[n], BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP,
             BHP[b].WindEnergyBuffer, eth, ratio, backlog, cap_hit);
      fflush(stdout);

    }

  bh_ffr_validate_state("post-wind-feedback");

  myfree(FeedbackWinnerID);
  myfree(FeedbackResults);
  myfree(FeedbackEvents);
  myfree(FeedbackPacketCount);

  FeedbackWinnerID = NULL;
  FeedbackResults = NULL;
  FeedbackEvents = NULL;
  FeedbackPacketCount = NULL;
  FeedbackNTargets = 0;
}
