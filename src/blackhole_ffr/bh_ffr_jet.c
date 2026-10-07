#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Iteration 8: persistent jet-axis memory and narrow bipolar jet feedback.
 *
 * JetDir is a numerical memory variable, not the physical BH spin. When the
 * unresolved reservoir is coherent, the jet/disc misalignment angle decays
 * exactly as theta(t+dt)=theta(t) exp[-dt/(xi_dir tau_d)]. When coherence is
 * below BHMinCoherence, JetDir is frozen.
 *
 * Jet feedback uses the same burst-threshold philosophy as the validated wind
 * channel, but carries no independent rest-mass sink or mass return. Active
 * gas cells in the two narrow JetDir cones receive separately normalized,
 * mass-weighted bipolar impulses. The impulse amplitude is solved from the
 * exact kinetic-energy quadratic in physical peculiar velocity, giving zero
 * net kick momentum and the requested packet energy to MPI roundoff.
 */

enum bh_ffr_jet_pass
{
  BH_FFR_JET_STATS = 0,
  BH_FFR_JET_MARK = 1,
  BH_FFR_JET_CONFLICT = 2,
  BH_FFR_JET_INJECT = 3
};

struct bh_ffr_jet_event
{
  int Candidate;
  int Fire;
  MyDouble PacketEnergy;
  MyDouble Q;
  MyDouble LobeMass[2];
  long long LobeCount[2];
  MyDouble EnergyBefore;
  MyFloat CosCone;
  int UsedHemisphereFallback;
};

static struct bh_ffr_jet_event *JetEvents;
static int *JetPacketCount;
static MyIDType *JetWinnerID;
static int JetNTargets;
static int JetPass;
static int JetConflictRound;

typedef struct
{
  MyDouble Pos[3];
  MyFloat Axis[3];
  MyFloat Radius;
  MyFloat LiveRadius;
  MyDouble AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT];
  MyFloat CosCone;
  MyIDType BHID;
  int Candidate;
  int Fire;
  MyDouble Q;
  MyDouble LobeMass[2];
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble EnclosedMass;
  MyDouble LobeMass[2];
  MyDouble LobeProjectedMomentum[2];
  MyDouble LobeTotalMass[2];
  long long LobeCount[2];

  MyDouble HemiMass[2];
  MyDouble HemiProjectedMomentum[2];
  MyDouble HemiTotalMass[2];
  long long HemiCount[2];

  /* Diagnostics-only adaptive MACER jet-radius survey. */
  MyDouble AdaptiveLobeMass[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  MyDouble AdaptiveLobeTotalMass[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  long long AdaptiveLobeCount[BH_FFR_ADAPTIVE_RADIUS_COUNT][2];
  MyDouble AdaptiveHemiMass[2];
  MyDouble AdaptiveHemiTotalMass[2];
  long long AdaptiveHemiCount[2];

  int Conflict;

  MyDouble KickEnergy;
  MyDouble KickMomentum[3];
} data_out;

static data_out *JetResults;
static data_out *DataResult, *DataOut;

static double bh_ffr_jet_coordinate_radius(void)
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;

  if(!isfinite(a) || !(a > 0) || !isfinite(All.BHFeedbackRadius) || !(All.BHFeedbackRadius > 0))
    terminate("BH_FFR: invalid jet feedback-radius conversion Rfb=%g a=%g", All.BHFeedbackRadius, a);

  const double r = All.BHFeedbackRadius / a;
  if(!isfinite(r) || !(r > 0))
    terminate("BH_FFR: invalid coordinate jet feedback radius=%g", r);

  return r;
}

static int bh_ffr_jet_decode_hydro_timebin(int bin)
{
  if(bin < 0)
    bin = -bin - 1;
  return bin;
}

static int bh_ffr_jet_gas_is_active(int j)
{
  const int bin = bh_ffr_jet_decode_hydro_timebin(P[j].TimeBinHydro);
  return bin > 0 && bin < TIMEBINS && TimeBinSynchronized[bin];
}

static int bh_ffr_jet_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= JetNTargets)
    terminate("BH_FFR: jet target=%d outside [0,%d) in %s", target, JetNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active jet particle=%d in %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact jet state for particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: jet particle ID=%llu is not drifted to Ti_Current=%lld in %s", (unsigned long long)P[p].ID,
              (long long)All.Ti_Current, where);

  return p;
}

static double bh_ffr_dot3(const MyDouble a[3], const MyDouble b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double bh_ffr_norm3(const MyDouble a[3])
{
  return sqrt(bh_ffr_dot3(a, a));
}

static void bh_ffr_normalize3(MyDouble a[3], const char *where)
{
  const double n = bh_ffr_norm3(a);
  if(!isfinite(n) || !(n > 0))
    terminate("BH_FFR: cannot normalize jet vector in %s", where);

  for(int k = 0; k < 3; k++)
    a[k] /= n;
}

/* Exact spherical interpolation for a fraction f of the great-circle angle.
 * The antiparallel case has no unique great circle, so choose a deterministic
 * perpendicular axis and rotate along that circle. */
static void bh_ffr_slerp_axis(const MyDouble from_in[3], const MyDouble to_in[3], double f, MyDouble out[3])
{
  if(!isfinite(f) || f < 0 || f > 1)
    terminate("BH_FFR: invalid jet slerp fraction=%g", f);

  MyDouble from[3] = {from_in[0], from_in[1], from_in[2]};
  MyDouble to[3] = {to_in[0], to_in[1], to_in[2]};
  bh_ffr_normalize3(from, "slerp/from");
  bh_ffr_normalize3(to, "slerp/to");

  double dot = bh_ffr_dot3(from, to);
  dot = dmax(-1.0, dmin(1.0, dot));

  if(f == 0)
    {
      for(int k = 0; k < 3; k++)
        out[k] = from[k];
      return;
    }
  if(f == 1)
    {
      for(int k = 0; k < 3; k++)
        out[k] = to[k];
      return;
    }

  if(dot > 1.0 - 1.0e-12)
    {
      for(int k = 0; k < 3; k++)
        out[k] = (1.0 - f) * from[k] + f * to[k];
      bh_ffr_normalize3(out, "slerp/parallel");
      return;
    }

  if(dot < -1.0 + 1.0e-10)
    {
      MyDouble ref[3] = {0, 0, 0};
      int imin = 0;
      if(fabs(from[1]) < fabs(from[imin]))
        imin = 1;
      if(fabs(from[2]) < fabs(from[imin]))
        imin = 2;
      ref[imin] = 1.0;

      MyDouble perp[3] = {
          from[1] * ref[2] - from[2] * ref[1],
          from[2] * ref[0] - from[0] * ref[2],
          from[0] * ref[1] - from[1] * ref[0]};
      bh_ffr_normalize3(perp, "slerp/antiparallel");

      const double angle = M_PI * f;
      for(int k = 0; k < 3; k++)
        out[k] = cos(angle) * from[k] + sin(angle) * perp[k];
      bh_ffr_normalize3(out, "slerp/antiparallel-result");
      return;
    }

  const double theta = acos(dot);
  const double sintheta = sin(theta);
  const double w0 = sin((1.0 - f) * theta) / sintheta;
  const double w1 = sin(f * theta) / sintheta;

  for(int k = 0; k < 3; k++)
    out[k] = w0 * from[k] + w1 * to[k];

  bh_ffr_normalize3(out, "slerp/general");
}

void bh_ffr_update_jet_direction(int p, double dt_myr, double disk_time_myr)
{
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid particle index=%d in jet-direction update", p);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for particle ID=%llu in jet-direction update", (unsigned long long)P[p].ID);

  if(!isfinite(dt_myr) || dt_myr < 0)
    terminate("BH_FFR: invalid jet-direction dt=%g Myr", dt_myr);

  double coherence_norm2 = 0.0;
  for(int k = 0; k < 3; k++)
    coherence_norm2 += BHP[b].Coherence[k] * BHP[b].Coherence[k];

  const double coherence =
      BHP[b].ReservoirMass > 0 ? sqrt(coherence_norm2) / BHP[b].ReservoirMass : 0.0;

  if(!isfinite(coherence) || coherence < 0 || coherence > 1.0 + 1.0e-10)
    terminate("BH_FFR: invalid jet-direction coherence=%g for ID=%llu", coherence, (unsigned long long)P[p].ID);

  if(dt_myr == 0 || coherence < All.BHMinCoherence)
    return;

  if(!isfinite(disk_time_myr) || !(disk_time_myr > 0))
    terminate("BH_FFR: invalid frozen reservoir time=%g Myr in jet-direction update", disk_time_myr);
  double tdir_myr = DBL_MAX;
  if(disk_time_myr < DBL_MAX / All.BHJetDirectionTimeFactor)
    tdir_myr = All.BHJetDirectionTimeFactor * disk_time_myr;
  if(!isfinite(tdir_myr) || !(tdir_myr > 0))
    terminate("BH_FFR: invalid jet alignment time=%g Myr", tdir_myr);

  MyDouble jet_old[3] = {BHP[b].JetDir[0], BHP[b].JetDir[1], BHP[b].JetDir[2]};
  MyDouble disc[3] = {BHP[b].DiscDir[0], BHP[b].DiscDir[1], BHP[b].DiscDir[2]};
  bh_ffr_normalize3(jet_old, "jet-update/old");
  bh_ffr_normalize3(disc, "jet-update/disc");

  double dot0 = bh_ffr_dot3(jet_old, disc);
  dot0 = dmax(-1.0, dmin(1.0, dot0));
  const double theta0 = acos(dot0);

  if(theta0 <= 1.0e-14)
    {
      for(int k = 0; k < 3; k++)
        BHP[b].JetDir[k] = disc[k];
      return;
    }

  const double decay = exp(-dt_myr / tdir_myr);
  const double f = 1.0 - decay;
  MyDouble updated[3];
  bh_ffr_slerp_axis(jet_old, disc, f, updated);

  double dot1 = bh_ffr_dot3(updated, disc);
  dot1 = dmax(-1.0, dmin(1.0, dot1));
  const double theta1 = acos(dot1);
  const double expected = theta0 * decay;

  if(fabs(theta1 - expected) > 2.0e-10 * fmax(theta0, 1.0e-30))
    terminate("BH_FFR: jet-axis exponential update failed ID=%llu theta1=%g expected=%g",
              (unsigned long long)P[p].ID, theta1, expected);

  for(int k = 0; k < 3; k++)
    BHP[b].JetDir[k] = updated[k];

  printf("BH_FFR: jet axis update ID=%llu task=%d theta0=%.17g theta1=%.17g dtMyr=%.17g tdirMyr=%.17g coherence=%.17g\n",
         (unsigned long long)P[p].ID, ThisTask, theta0, theta1, dt_myr, tdir_myr, coherence);
  fflush(stdout);
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_jet_particle_from_target(target, "particle2in");
  const int b = P[p].BHDataIndex;
  const struct bh_ffr_jet_event *ev = &JetEvents[target];

  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[p].Pos[k];
      in->Axis[k] = BHP[b].JetDir[k];
    }

  in->LiveRadius = bh_ffr_jet_coordinate_radius();
  in->Radius = in->LiveRadius;
  for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
    in->AdaptiveRadius[k] = 0.0;

  if(JetPass == BH_FFR_JET_STATS &&
     All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_MACER)
    {
      struct bh_ffr_discrete_radius_grid grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      if(grid.Count != BH_FFR_ADAPTIVE_RADIUS_COUNT)
        terminate("BH_FFR: unexpected jet adaptive-radius count=%d", grid.Count);

      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0))
        terminate("BH_FFR: invalid scale factor=%g in jet adaptive-radius setup", a);

      for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
        in->AdaptiveRadius[k] = grid.Radius[k] / a;

      if(in->AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT - 1] > in->Radius)
        in->Radius = in->AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT - 1];
    }

  in->CosCone = (JetPass == BH_FFR_JET_STATS)
                    ? cos(All.BHJetConeAngleDeg * M_PI / 180.0)
                    : ev->CosCone;
  in->BHID = P[p].ID;
  in->Candidate = ev->Candidate;
  in->Fire = ev->Fire;
  in->Q = ev->Q;
  in->LobeMass[0] = ev->LobeMass[0];
  in->LobeMass[1] = ev->LobeMass[1];
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= JetNTargets)
    terminate("BH_FFR: jet result target=%d outside [0,%d)", target, JetNTargets);

  data_out *res = &JetResults[target];

  if(mode == MODE_LOCAL_PARTICLES)
    *res = *out;
  else
    {
      res->EnclosedMass += out->EnclosedMass;
      for(int l = 0; l < 2; l++)
        {
          res->LobeMass[l] += out->LobeMass[l];
          res->LobeProjectedMomentum[l] += out->LobeProjectedMomentum[l];
          res->LobeTotalMass[l] += out->LobeTotalMass[l];
          res->LobeCount[l] += out->LobeCount[l];

          res->HemiMass[l] += out->HemiMass[l];
          res->HemiProjectedMomentum[l] += out->HemiProjectedMomentum[l];
          res->HemiTotalMass[l] += out->HemiTotalMass[l];
          res->HemiCount[l] += out->HemiCount[l];

          res->AdaptiveHemiMass[l] += out->AdaptiveHemiMass[l];
          res->AdaptiveHemiTotalMass[l] += out->AdaptiveHemiTotalMass[l];
          res->AdaptiveHemiCount[l] += out->AdaptiveHemiCount[l];
        }
      for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
        for(int l = 0; l < 2; l++)
          {
            res->AdaptiveLobeMass[k][l] += out->AdaptiveLobeMass[k][l];
            res->AdaptiveLobeTotalMass[k][l] += out->AdaptiveLobeTotalMass[k][l];
            res->AdaptiveLobeCount[k][l] += out->AdaptiveLobeCount[k][l];
          }

      if(out->Conflict)
        res->Conflict = 1;

      res->KickEnergy += out->KickEnergy;
      for(int k = 0; k < 3; k++)
        res->KickMomentum[k] += out->KickMomentum[k];
    }
}

#include "../utils/generic_comm_helpers2.h"

static int bh_ffr_jet_evaluate(int target, int mode, int threadid);

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
      if(target >= JetNTargets)
        break;

      bh_ffr_jet_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_jet_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static unsigned long long bh_ffr_jet_priority(MyIDType id)
{
  unsigned long long x = (unsigned long long)id;
  x ^= (unsigned long long)All.Ti_Current + 0x9e3779b97f4a7c15ULL +
       (unsigned long long)(JetConflictRound + 1) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static int bh_ffr_jet_candidate_wins(MyIDType candidate, MyIDType incumbent)
{
  if(incumbent == 0)
    return 1;

  const unsigned long long pc = bh_ffr_jet_priority(candidate);
  const unsigned long long pi = bh_ffr_jet_priority(incumbent);
  return pc < pi || (pc == pi && candidate < incumbent);
}

static int bh_ffr_jet_selected_lobe(const data_in *in, int j, double *r2_out)
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

static int bh_ffr_jet_evaluate(int target, int mode, int threadid)
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

  if(JetPass != BH_FFR_JET_STATS)
    {
      if((JetPass == BH_FFR_JET_INJECT && !in->Fire) ||
         ((JetPass == BH_FFR_JET_MARK || JetPass == BH_FFR_JET_CONFLICT) && !in->Candidate))
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
    terminate("BH_FFR: invalid scale factor=%g in jet evaluation", a);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_FFR: jet gas index=%d outside NumGas=%d", j, NumGas);
      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      double r2;
      const int lobe = bh_ffr_jet_selected_lobe(in, j, &r2);

      if(JetPass == BH_FFR_JET_STATS)
        {
          const int active = bh_ffr_jet_gas_is_active(j);
          int hemi_search = -1;

          if(r2 > 0 && r2 <= in->Radius * in->Radius)
            {
              const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
              const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
              const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
              const double dot = dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
              hemi_search = dot >= 0 ? 0 : 1;
            }

          if(hemi_search >= 0 &&
             in->AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT - 1] > 0)
            {
              for(int k = 0; k < BH_FFR_ADAPTIVE_RADIUS_COUNT; k++)
                if(r2 <= in->AdaptiveRadius[k] * in->AdaptiveRadius[k] && lobe >= 0)
                  {
                    out.AdaptiveLobeTotalMass[k][lobe] += P[j].Mass;
                    if(active)
                      {
                        out.AdaptiveLobeCount[k][lobe]++;
                        out.AdaptiveLobeMass[k][lobe] += P[j].Mass;
                      }
                  }

              const double rmax =
                  in->AdaptiveRadius[BH_FFR_ADAPTIVE_RADIUS_COUNT - 1];
              if(r2 <= rmax * rmax)
                {
                  out.AdaptiveHemiTotalMass[hemi_search] += P[j].Mass;
                  if(active)
                    {
                      out.AdaptiveHemiCount[hemi_search]++;
                      out.AdaptiveHemiMass[hemi_search] += P[j].Mass;
                    }
                }
            }

          int hemi = -1;
          if(r2 > 0 && r2 <= in->LiveRadius * in->LiveRadius)
            {
              hemi = hemi_search;
              out.EnclosedMass += P[j].Mass;
              out.HemiTotalMass[hemi] += P[j].Mass;
              if(lobe >= 0)
                out.LobeTotalMass[lobe] += P[j].Mass;
            }

          if(active && hemi >= 0)
            {
              double vdot = 0.0;
              for(int k = 0; k < 3; k++)
                vdot += (P[j].Vel[k] / a) * in->Axis[k];

              out.HemiCount[hemi]++;
              out.HemiMass[hemi] += P[j].Mass;
              out.HemiProjectedMomentum[hemi] += P[j].Mass * vdot;

              if(lobe >= 0)
                {
                  out.LobeCount[lobe]++;
                  out.LobeMass[lobe] += P[j].Mass;
                  out.LobeProjectedMomentum[lobe] += P[j].Mass * vdot;
                }
            }
          continue;
        }

      if(lobe < 0 || !bh_ffr_jet_gas_is_active(j))
        continue;

      if(JetPass == BH_FFR_JET_MARK)
        {
          if(bh_ffr_jet_candidate_wins(in->BHID, JetWinnerID[j]))
            JetWinnerID[j] = in->BHID;
          continue;
        }

      if(JetPass == BH_FFR_JET_CONFLICT)
        {
          if(JetWinnerID[j] != in->BHID)
            out.Conflict = 1;
          continue;
        }

      if(JetPass != BH_FFR_JET_INJECT)
        terminate("BH_FFR: unknown jet pass=%d", JetPass);

      if(!(in->LobeMass[lobe] > 0))
        terminate("BH_FFR: non-positive jet lobe mass for BH ID=%llu", (unsigned long long)in->BHID);

      const double weight = P[j].Mass / in->LobeMass[lobe];
      if(!isfinite(weight) || !(weight > 0))
        terminate("BH_FFR: invalid jet target weight=%g for gas ID=%llu", weight, (unsigned long long)P[j].ID);

      const double sign = (lobe == 0) ? 1.0 : -1.0;
      double dp_phys[3], v_phys[3];
      double dp2 = 0.0, vdotdp = 0.0;

      for(int k = 0; k < 3; k++)
        {
          v_phys[k] = P[j].Vel[k] / a;
          dp_phys[k] = sign * in->Q * weight * in->Axis[k];
          dp2 += dp_phys[k] * dp_phys[k];
          vdotdp += v_phys[k] * dp_phys[k];
        }

      const double dkin = vdotdp + 0.5 * dp2 / P[j].Mass;
      if(!isfinite(dkin))
        terminate("BH_FFR: non-finite jet kinetic increment for gas ID=%llu", (unsigned long long)P[j].ID);

      for(int k = 0; k < 3; k++)
        {
          SphP[j].Momentum[k] += a * dp_phys[k];
          out.KickMomentum[k] += dp_phys[k];
        }

      SphP[j].Energy += a * a * dkin;

      for(int k = 0; k < 3; k++)
        P[j].Vel[k] = SphP[j].Momentum[k] / P[j].Mass;

      out.KickEnergy += dkin;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static void bh_ffr_jet_comm_pass(int pass)
{
  JetPass = pass;
  memset(JetResults, 0, (JetNTargets > 0 ? JetNTargets : 1) * sizeof(*JetResults));

  if(All.TotNumGas > 0)
    {
      generic_set_MaxNexport();
      generic_comm_pattern(JetNTargets, kernel_local, kernel_imported);
    }
}

static double bh_ffr_jet_packet_q(double A, double B, double energy)
{
  if(!isfinite(A) || !isfinite(B) || !(B > 0) || !isfinite(energy) || !(energy > 0))
    terminate("BH_FFR: invalid jet packet-root inputs A=%g B=%g E=%g", A, B, energy);

  const double disc = sqrt(A * A + 2.0 * B * energy);
  const double q = (A >= 0) ? 2.0 * energy / (disc + A) : (-A + disc) / B;

  if(!isfinite(q) || !(q > 0))
    terminate("BH_FFR: invalid positive jet packet root q=%g A=%g B=%g E=%g", q, A, B, energy);

  return q;
}

static int bh_ffr_jet_survey_radius_ok(const data_out *res, int k,
                                        double *fplus, double *fminus)
{
  const double mtot_plus = res->AdaptiveLobeTotalMass[k][0];
  const double mtot_minus = res->AdaptiveLobeTotalMass[k][1];
  *fplus = mtot_plus > 0 ? res->AdaptiveLobeMass[k][0] / mtot_plus : 0.0;
  *fminus = mtot_minus > 0 ? res->AdaptiveLobeMass[k][1] / mtot_minus : 0.0;

  return res->AdaptiveLobeCount[k][0] >= 4 &&
         res->AdaptiveLobeCount[k][1] >= 4 &&
         mtot_plus > 0 && mtot_minus > 0 &&
         *fplus >= All.BHMinActiveTargetMassFrac &&
         *fminus >= All.BHMinActiveTargetMassFrac;
}

static void bh_ffr_jet_report_adaptive_radius_survey(void)
{
  if(All.BHBenchmarkFeedbackModel != BH_BENCHMARK_FEEDBACK_MACER ||
     JetConflictRound != 0)
    return;

  static unsigned long long survey_calls = 0;
  if((survey_calls++ % 64ULL) != 0ULL)
    return;

  for(int n = 0; n < JetNTargets; n++)
    {
      const int p =
          bh_ffr_jet_particle_from_target(n, "bh_ffr_jet_report_adaptive_radius_survey");
      const data_out *res = &JetResults[n];

      struct bh_ffr_discrete_radius_grid grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      const double eps = bh_ffr_effective_softening_proper(p);

      int selected = -1;
      double fplus = 0.0, fminus = 0.0;
      for(int k = 0; k < grid.Count; k++)
        {
          double fp = 0.0, fm = 0.0;
          const int ok = bh_ffr_jet_survey_radius_ok(res, k, &fp, &fm);
          printf("BH_FFR: adaptive jet radius detail mode=survey ID=%llu task=%d "
                 "index=%d R=%g Reps=%g Nplus=%lld Nminus=%lld "
                 "fplus=%g fminus=%g resolved=%d\n",
                 (unsigned long long)P[p].ID, ThisTask, k,
                 (double)grid.Radius[k], (double)(grid.Radius[k] / eps),
                 res->AdaptiveLobeCount[k][0], res->AdaptiveLobeCount[k][1],
                 fp, fm, ok);
          if(selected < 0 && ok)
            {
              selected = k;
              fplus = fp;
              fminus = fm;
            }
        }

      int fallback = 0;
      int underresolved = 0;
      if(selected < 0)
        {
          const double mtot_plus = res->AdaptiveHemiTotalMass[0];
          const double mtot_minus = res->AdaptiveHemiTotalMass[1];
          fplus = mtot_plus > 0 ? res->AdaptiveHemiMass[0] / mtot_plus : 0.0;
          fminus = mtot_minus > 0 ? res->AdaptiveHemiMass[1] / mtot_minus : 0.0;
          const int hemi_ok =
              res->AdaptiveHemiCount[0] >= 4 &&
              res->AdaptiveHemiCount[1] >= 4 &&
              mtot_plus > 0 && mtot_minus > 0 &&
              fplus >= All.BHMinActiveTargetMassFrac &&
              fminus >= All.BHMinActiveTargetMassFrac;
          selected = grid.Count - 1;
          fallback = hemi_ok ? 1 : 0;
          underresolved = hemi_ok ? 0 : 1;
        }

      const long long nplus = fallback ? res->AdaptiveHemiCount[0]
                                       : res->AdaptiveLobeCount[selected][0];
      const long long nminus = fallback ? res->AdaptiveHemiCount[1]
                                        : res->AdaptiveLobeCount[selected][1];

      printf("BH_FFR: adaptive jet radius mode=survey ID=%llu task=%d eps=%g "
             "Rsearch=%g Rjet=%g index=%d Nplus=%lld Nminus=%lld "
             "fplus=%g fminus=%g fallback=%d underresolved=%d\n",
             (unsigned long long)P[p].ID, ThisTask, eps,
             (double)grid.Radius[grid.Count - 1],
             (double)grid.Radius[selected], selected, nplus, nminus,
             fplus, fminus, fallback, underresolved);
      fflush(stdout);
    }
}

static void bh_ffr_jet_prepare_candidates(void)
{
  for(int n = 0; n < JetNTargets; n++)
    {
      struct bh_ffr_jet_event *ev = &JetEvents[n];
      memset(ev, 0, sizeof(*ev));

      const int p = bh_ffr_jet_particle_from_target(n, "bh_ffr_jet_prepare_candidates");
      const int b = P[p].BHDataIndex;
      const data_out *res = &JetResults[n];

      const double sigma2 = BHP[b].SigmaDM * BHP[b].SigmaDM;
      double vbind2 = sigma2;
      if(All.BHUseCentralBindingTerm)
        vbind2 += 2.0 * All.G * bh_ffr_central_mass_code(p) / All.BHFeedbackRadius;

      if(!isfinite(vbind2) || vbind2 < 0 || !isfinite(res->EnclosedMass) || res->EnclosedMass < 0)
        terminate("BH_FFR: invalid jet threshold environment Menc=%g vbind2=%g for ID=%llu", res->EnclosedMass, vbind2,
                  (unsigned long long)P[p].ID);

      const double base_threshold = 0.5 * res->EnclosedMass * vbind2;
      /* Keep wind/jet threshold ownership independent.  The live adaptive
       * wind pass computes WindThresholdEnergy from its own selected radius. */
      BHP[b].JetThresholdEnergy = All.BHJetBurstFactor * base_threshold;

      if(JetPacketCount[n] >= All.BHMaxPacketsPerStep)
        continue;

      const double eth = BHP[b].JetThresholdEnergy;
      if(!(eth > 0) || !(BHP[b].JetEnergyBuffer >= eth))
        continue;

      const double cone_active_mass = res->LobeMass[0] + res->LobeMass[1];
      const double cone_total_mass = res->LobeTotalMass[0] + res->LobeTotalMass[1];
      const int cone_ok =
          res->LobeCount[0] >= All.BHMinTargetsPerLobe &&
          res->LobeCount[1] >= All.BHMinTargetsPerLobe &&
          res->LobeMass[0] > 0 && res->LobeMass[1] > 0 &&
          cone_total_mass > 0 &&
          cone_active_mass >= All.BHMinActiveTargetMassFrac * cone_total_mass;

      const double hemi_active_mass = res->HemiMass[0] + res->HemiMass[1];
      const double hemi_total_mass = res->HemiTotalMass[0] + res->HemiTotalMass[1];
      const int hemi_ok =
          res->HemiCount[0] >= All.BHMinTargetsPerLobe &&
          res->HemiCount[1] >= All.BHMinTargetsPerLobe &&
          res->HemiMass[0] > 0 && res->HemiMass[1] > 0 &&
          hemi_total_mass > 0 &&
          hemi_active_mass >= All.BHMinActiveTargetMassFrac * hemi_total_mass;

      if(!cone_ok && !hemi_ok)
        continue;

      ev->UsedHemisphereFallback = cone_ok ? 0 : 1;
      ev->CosCone = cone_ok ? cos(All.BHJetConeAngleDeg * M_PI / 180.0) : 0.0;

      ev->EnergyBefore = BHP[b].JetEnergyBuffer;

      /* Release one binding-threshold quantum per packet.  Excess jet
       * energy remains buffered for later rounds/steps; this controls burst
       * strength without reintroducing feedback-driven BH subcycling. */
      ev->PacketEnergy = eth;
      for(int l = 0; l < 2; l++)
        {
          ev->LobeMass[l] = ev->UsedHemisphereFallback ? res->HemiMass[l] : res->LobeMass[l];
          ev->LobeCount[l] = ev->UsedHemisphereFallback ? res->HemiCount[l] : res->LobeCount[l];
        }

      const double pplus =
          ev->UsedHemisphereFallback ? res->HemiProjectedMomentum[0] : res->LobeProjectedMomentum[0];
      const double pminus =
          ev->UsedHemisphereFallback ? res->HemiProjectedMomentum[1] : res->LobeProjectedMomentum[1];
      if(!isfinite(ev->LobeMass[0]) || !isfinite(ev->LobeMass[1]) ||
         !(ev->LobeMass[0] > 0) || !(ev->LobeMass[1] > 0))
        terminate("BH_FFR: invalid selected jet lobe masses ID=%llu fallback=%d Mplus=%g Mminus=%g",
                  (unsigned long long)P[p].ID, ev->UsedHemisphereFallback,
                  (double)ev->LobeMass[0], (double)ev->LobeMass[1]);

      const double A =
          pplus / ev->LobeMass[0] - pminus / ev->LobeMass[1];
      const double B = 1.0 / ev->LobeMass[0] + 1.0 / ev->LobeMass[1];

      ev->Q = bh_ffr_jet_packet_q(A, B, ev->PacketEnergy);
      ev->Candidate = 1;
    }
}

static int bh_ffr_jet_global_candidate_count(void)
{
  int local = 0;
  for(int n = 0; n < JetNTargets; n++)
    if(JetEvents[n].Candidate)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static int bh_ffr_jet_global_fire_count(void)
{
  int local = 0;
  for(int n = 0; n < JetNTargets; n++)
    if(JetEvents[n].Fire)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static void bh_ffr_jet_resolve_conflicts(void)
{
  const int candidates = bh_ffr_jet_global_candidate_count();
  if(candidates <= 0)
    return;

  if(candidates == 1)
    {
      for(int n = 0; n < JetNTargets; n++)
        JetEvents[n].Fire = JetEvents[n].Candidate;
      return;
    }

  memset(JetWinnerID, 0, (NumGas > 0 ? NumGas : 1) * sizeof(*JetWinnerID));
  bh_ffr_jet_comm_pass(BH_FFR_JET_MARK);
  bh_ffr_jet_comm_pass(BH_FFR_JET_CONFLICT);

  for(int n = 0; n < JetNTargets; n++)
    JetEvents[n].Fire = JetEvents[n].Candidate && !JetResults[n].Conflict;
}

static void bh_ffr_jet_commit_packets(void)
{
  for(int n = 0; n < JetNTargets; n++)
    {
      struct bh_ffr_jet_event *ev = &JetEvents[n];
      if(!ev->Fire)
        continue;

      const int p = bh_ffr_jet_particle_from_target(n, "bh_ffr_jet_commit_packets");
      const int b = P[p].BHDataIndex;
      const data_out *res = &JetResults[n];

      const double etol = 3.0e-8 * fmax(fabs(ev->PacketEnergy), 1.0e-30);
      if(fabs(res->KickEnergy - ev->PacketEnergy) > etol)
        terminate("BH_FFR: jet packet kinetic-energy mismatch ID=%llu injected=%g requested=%g",
                  (unsigned long long)P[p].ID, res->KickEnergy, ev->PacketEnergy);

      double pnorm2 = 0.0;
      for(int k = 0; k < 3; k++)
        pnorm2 += res->KickMomentum[k] * res->KickMomentum[k];
      const double pnorm = sqrt(pnorm2);

      if(pnorm > 3.0e-8 * fmax(fabs(ev->Q), 1.0e-30))
        terminate("BH_FFR: jet packet kick-momentum imbalance ID=%llu |sum dp|=%g q=%g", (unsigned long long)P[p].ID,
                  pnorm, ev->Q);

      BHP[b].JetEnergyBuffer = ev->EnergyBefore - ev->PacketEnergy;
      if(BHP[b].JetEnergyBuffer < 0 && BHP[b].JetEnergyBuffer > -3.0e-12 * fmax(fabs(ev->EnergyBefore), 1.0e-30))
        BHP[b].JetEnergyBuffer = 0.0;
      if(BHP[b].JetEnergyBuffer < 0)
        terminate("BH_FFR: negative jet buffer after packet for ID=%llu", (unsigned long long)P[p].ID);

      JetPacketCount[n]++;

      const double eerg = ev->PacketEnergy * All.UnitEnergy_in_cgs / All.HubbleParam;
      printf("BH_FFR: jet packet fired ID=%llu task=%d E=%g erg q=%g Nplus=%lld Nminus=%lld dErel=%g pbal=%g\n",
             (unsigned long long)P[p].ID, ThisTask, eerg, ev->Q, ev->LobeCount[0], ev->LobeCount[1],
             fabs(res->KickEnergy - ev->PacketEnergy) / ev->PacketEnergy, pnorm / fmax(fabs(ev->Q), 1.0e-30));
      printf("BH_FFR: jet coupling ID=%llu task=%d fallback=%d cosCone=%g Nplus=%lld Nminus=%lld\n",
             (unsigned long long)P[p].ID, ThisTask, ev->UsedHemisphereFallback,
             (double)ev->CosCone, ev->LobeCount[0], ev->LobeCount[1]);
      fflush(stdout);
    }
}

void bh_ffr_jet_self_test(void)
{
  MyDouble from[3] = {1.0, 0.0, 0.0};
  MyDouble to[3] = {0.0, 1.0, 0.0};
  MyDouble out[3];

  const double decay = exp(-1.0);
  bh_ffr_slerp_axis(from, to, 1.0 - decay, out);

  double dot = bh_ffr_dot3(out, to);
  dot = dmax(-1.0, dmin(1.0, dot));
  const double theta = acos(dot);
  const double expected_theta = 0.5 * M_PI * decay;

  if(fabs(theta - expected_theta) > 2.0e-13)
    terminate("BH_FFR: jet self-test failed exponential slerp theta=%g expected=%g", theta, expected_theta);

  const double A = 0.23;
  const double B = 0.77;
  const double energy = 1.17;
  const double q = bh_ffr_jet_packet_q(A, B, energy);
  const double got = A * q + 0.5 * B * q * q;

  if(fabs(got - energy) > 2.0e-13 * energy)
    terminate("BH_FFR: jet self-test failed packet root got=%g expected=%g", got, energy);
}

void bh_ffr_inject_jet_feedback(void)
{
  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

#ifdef MHD
  terminate("BH_FFR: kinetic jet feedback is not enabled with MHD until the magnetic feedback policy is validated");
#endif

  JetNTargets = NumActiveBHFFR;
  JetPacketCount = (int *)mymalloc("BHFFRJetPacketCount", (JetNTargets > 0 ? JetNTargets : 1) * sizeof(int));
  JetEvents =
      (struct bh_ffr_jet_event *)mymalloc("BHFFRJetEvents", (JetNTargets > 0 ? JetNTargets : 1) * sizeof(*JetEvents));
  JetResults = (data_out *)mymalloc("BHFFRJetResults", (JetNTargets > 0 ? JetNTargets : 1) * sizeof(*JetResults));
  JetWinnerID = (MyIDType *)mymalloc("BHFFRJetWinner", (NumGas > 0 ? NumGas : 1) * sizeof(*JetWinnerID));

  memset(JetPacketCount, 0, (JetNTargets > 0 ? JetNTargets : 1) * sizeof(int));
  memset(JetEvents, 0, (JetNTargets > 0 ? JetNTargets : 1) * sizeof(*JetEvents));

  int any_global_fire = 0;

  for(int round = 0; round < All.BHMaxPacketsPerStep; round++)
    {
      JetConflictRound = round;
      bh_ffr_jet_comm_pass(BH_FFR_JET_STATS);
      bh_ffr_jet_report_adaptive_radius_survey();
      bh_ffr_jet_prepare_candidates();

      if(bh_ffr_jet_global_candidate_count() <= 0)
        break;

      bh_ffr_jet_resolve_conflicts();
      const int global_fire = bh_ffr_jet_global_fire_count();
      if(global_fire <= 0)
        break;

      bh_ffr_jet_comm_pass(BH_FFR_JET_INJECT);
      bh_ffr_jet_commit_packets();
      any_global_fire += global_fire;
    }

  if(any_global_fire > 0)
    update_primitive_variables();

  bh_ffr_validate_state("post-jet-feedback");

  myfree(JetWinnerID);
  myfree(JetResults);
  myfree(JetEvents);
  myfree(JetPacketCount);

  JetWinnerID = NULL;
  JetResults = NULL;
  JetEvents = NULL;
  JetPacketCount = NULL;
  JetNTargets = 0;
}
