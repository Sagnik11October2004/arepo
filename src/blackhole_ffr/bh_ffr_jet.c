#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

#define BH_FFR_JET_BROADEN_ANGLE_COUNT 4

static const double BHFFRJetBroadenAnglesDeg[BH_FFR_JET_BROADEN_ANGLE_COUNT] =
    {15.0, 30.0, 45.0, 60.0};

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
  BH_FFR_JET_INJECT = 3,
  BH_FFR_JET_WAKE = 4
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
  MyDouble RadiusProper;
  MyFloat RadiusCoordinate;
  MyFloat CosCone;
  MyFloat AngleDeg;
  int AngleIndex;
  int Underresolved;
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
  MyDouble BroadenRadius[BH_FFR_JET_BROADEN_ANGLE_COUNT];
  MyFloat BroadenCosCone[BH_FFR_JET_BROADEN_ANGLE_COUNT];
  MyFloat CosCone;
  MyIDType BHID;
  int Candidate;
  int Fire;
  MyDouble Q;
  MyDouble LobeMass[2];
  int WakeTimeBin;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  /* Constant-volume jet broadening statistics. */
  MyDouble BroadenEnclosedMass[BH_FFR_JET_BROADEN_ANGLE_COUNT];
  MyDouble BroadenLobeMass[BH_FFR_JET_BROADEN_ANGLE_COUNT][2];
  MyDouble BroadenLobeProjectedMomentum[BH_FFR_JET_BROADEN_ANGLE_COUNT][2];
  MyDouble BroadenLobeTotalMass[BH_FFR_JET_BROADEN_ANGLE_COUNT][2];
  long long BroadenLobeCount[BH_FFR_JET_BROADEN_ANGLE_COUNT][2];

  int Conflict;

  MyDouble KickEnergy;
  MyDouble KickMomentum[3];
} data_out;

static data_out *JetResults;
static data_out *DataResult, *DataOut;

static int bh_ffr_jet_gas_is_active(int j)
{
  return bh_ffr_feedback_gas_is_hydro_active_now(j);
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

  in->Radius = 0.0;
  for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
    {
      in->BroadenRadius[q] = 0.0;
      in->BroadenCosCone[q] = 0.0;
    }

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

      const double theta_ref =
          BHFFRJetBroadenAnglesDeg[0] * M_PI / 180.0;
      const double omega_ref = 1.0 - cos(theta_ref);
      const double rref = grid.Radius[grid.Count - 1];

      for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
        {
          const double theta =
              BHFFRJetBroadenAnglesDeg[q] * M_PI / 180.0;
          const double omega = 1.0 - cos(theta);
          if(!isfinite(omega) || !(omega > 0))
            terminate("BH_FFR: invalid broadened jet solid-angle factor q=%d theta=%g",
                      q, BHFFRJetBroadenAnglesDeg[q]);

          const double rproper =
              rref * cbrt(omega_ref / omega);
          if(!isfinite(rproper) || !(rproper > 0))
            terminate("BH_FFR: invalid broadened jet radius q=%d R=%g",
                      q, rproper);

          in->BroadenRadius[q] = rproper / a;
          in->BroadenCosCone[q] = cos(theta);
        }

      in->Radius = in->BroadenRadius[0];
    }
  else
    {
      if(!isfinite(ev->RadiusCoordinate) || !(ev->RadiusCoordinate > 0))
        terminate("BH_FFR: invalid selected jet coordinate radius=%g for ID=%llu",
                  (double)ev->RadiusCoordinate, (unsigned long long)P[p].ID);
      in->Radius = ev->RadiusCoordinate;
    }

  in->CosCone = (JetPass == BH_FFR_JET_STATS)
                    ? cos(BHFFRJetBroadenAnglesDeg[0] * M_PI / 180.0)
                    : ev->CosCone;
  in->BHID = P[p].ID;
  in->WakeTimeBin =
      JetPass == BH_FFR_JET_WAKE
          ? bh_ffr_feedback_jet_target_hydro_timebin(p)
          : -1;
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
      for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
        {
          res->BroadenEnclosedMass[q] += out->BroadenEnclosedMass[q];
          for(int l = 0; l < 2; l++)
            {
              res->BroadenLobeMass[q][l] += out->BroadenLobeMass[q][l];
              res->BroadenLobeProjectedMomentum[q][l] +=
                  out->BroadenLobeProjectedMomentum[q][l];
              res->BroadenLobeTotalMass[q][l] += out->BroadenLobeTotalMass[q][l];
              res->BroadenLobeCount[q][l] += out->BroadenLobeCount[q][l];
            }
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

      if(JetPass == BH_FFR_JET_WAKE)
        {
          if(lobe >= 0 && in->WakeTimeBin > 0)
            bh_ffr_feedback_request_wakeup(j, in->WakeTimeBin);
          continue;
        }

      if(JetPass == BH_FFR_JET_STATS)
        {
          const int active = bh_ffr_jet_gas_is_active(j);

          if(r2 > 0 && r2 <= in->Radius * in->Radius)
            {
              /* Constant-volume live geometry statistics. */
              const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
              const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
              const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
              const double dot =
                  dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
              const double cosine = fabs(dot) / sqrt(r2);
              const int broaden_lobe = dot >= 0 ? 0 : 1;

              double broaden_vdot = 0.0;
              if(active)
                for(int k = 0; k < 3; k++)
                  broaden_vdot += (P[j].Vel[k] / a) * in->Axis[k];

              for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
                if(in->BroadenRadius[q] > 0 &&
                   r2 <= in->BroadenRadius[q] * in->BroadenRadius[q])
                  {
                    out.BroadenEnclosedMass[q] += P[j].Mass;

                    if(cosine >= in->BroadenCosCone[q])
                      {
                        out.BroadenLobeTotalMass[q][broaden_lobe] += P[j].Mass;
                        if(active)
                          {
                            out.BroadenLobeCount[q][broaden_lobe]++;
                            out.BroadenLobeMass[q][broaden_lobe] += P[j].Mass;
                            out.BroadenLobeProjectedMomentum[q][broaden_lobe] +=
                                P[j].Mass * broaden_vdot;
                          }
                      }
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

static int bh_ffr_jet_broaden_angle_ok(const data_out *res, int q,
                                        double *fplus, double *fminus)
{
  const double mtot_plus = res->BroadenLobeTotalMass[q][0];
  const double mtot_minus = res->BroadenLobeTotalMass[q][1];
  *fplus = mtot_plus > 0 ? res->BroadenLobeMass[q][0] / mtot_plus : 0.0;
  *fminus = mtot_minus > 0 ? res->BroadenLobeMass[q][1] / mtot_minus : 0.0;

  return res->BroadenLobeCount[q][0] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         res->BroadenLobeCount[q][1] >= BH_FFR_FEEDBACK_MIN_ACTIVE_PER_LOBE &&
         mtot_plus > 0 && mtot_minus > 0 &&
         *fplus >= All.BHMinActiveTargetMassFrac &&
         *fminus >= All.BHMinActiveTargetMassFrac;
}

static int bh_ffr_jet_select_live_broadening(const data_out *res,
                                              double *fplus,
                                              double *fminus)
{
  *fplus = 0.0;
  *fminus = 0.0;

  for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
    {
      double fp = 0.0, fm = 0.0;
      if(bh_ffr_jet_broaden_angle_ok(res, q, &fp, &fm))
        {
          *fplus = fp;
          *fminus = fm;
          return q;
        }
    }

  return -1;
}

static void bh_ffr_jet_report_broadening_survey(void)
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
          bh_ffr_jet_particle_from_target(
              n, "bh_ffr_jet_report_broadening_survey");
      const data_out *res = &JetResults[n];

      struct bh_ffr_discrete_radius_grid grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      const double eps = bh_ffr_effective_softening_proper(p);
      const double theta_ref =
          BHFFRJetBroadenAnglesDeg[0] * M_PI / 180.0;
      const double omega_ref = 1.0 - cos(theta_ref);
      const double rref = grid.Radius[grid.Count - 1];

      double selected_fplus = 0.0, selected_fminus = 0.0;
      const int selected =
          bh_ffr_jet_select_live_broadening(
              res, &selected_fplus, &selected_fminus);

      for(int q = 0; q < BH_FFR_JET_BROADEN_ANGLE_COUNT; q++)
        {
          const double theta =
              BHFFRJetBroadenAnglesDeg[q] * M_PI / 180.0;
          const double rproper =
              rref * cbrt(omega_ref / (1.0 - cos(theta)));

          double fp = 0.0, fm = 0.0;
          const int ok =
              bh_ffr_jet_broaden_angle_ok(res, q, &fp, &fm);

          printf("BH_FFR: jet broadening detail mode=live ID=%llu task=%d "
                 "angle=%g R=%g Reps=%g volumeRatio=1 "
                 "Nplus=%lld Nminus=%lld fplus=%g fminus=%g resolved=%d\n",
                 (unsigned long long)P[p].ID, ThisTask,
                 BHFFRJetBroadenAnglesDeg[q], rproper, rproper / eps,
                 res->BroadenLobeCount[q][0],
                 res->BroadenLobeCount[q][1],
                 fp, fm, ok);

        }

      if(selected >= 0)
        {
          const double theta =
              BHFFRJetBroadenAnglesDeg[selected] * M_PI / 180.0;
          const double rproper =
              rref * cbrt(omega_ref / (1.0 - cos(theta)));
          printf("BH_FFR: jet broadening mode=live ID=%llu task=%d "
                 "selectedAngle=%g Rjet=%g Reps=%g angleIndex=%d "
                 "Nplus=%lld Nminus=%lld fplus=%g fminus=%g underresolved=0\n",
                 (unsigned long long)P[p].ID, ThisTask,
                 BHFFRJetBroadenAnglesDeg[selected], rproper, rproper / eps,
                 selected,
                 res->BroadenLobeCount[selected][0],
                 res->BroadenLobeCount[selected][1],
                 selected_fplus, selected_fminus);
        }
      else
        {
          printf("BH_FFR: jet broadening mode=live ID=%llu task=%d "
                 "selectedAngle=none Rjet=0 Reps=0 angleIndex=-1 "
                 "Nplus=0 Nminus=0 fplus=0 fminus=0 underresolved=1\n",
                 (unsigned long long)P[p].ID, ThisTask);
        }
      fflush(stdout);
    }
}

static void bh_ffr_jet_prepare_candidates(void)
{
  for(int n = 0; n < JetNTargets; n++)
    {
      struct bh_ffr_jet_event *ev = &JetEvents[n];
      memset(ev, 0, sizeof(*ev));

      const int p =
          bh_ffr_jet_particle_from_target(
              n, "bh_ffr_jet_prepare_candidates");
      const int b = P[p].BHDataIndex;
      const data_out *res = &JetResults[n];

      struct bh_ffr_discrete_radius_grid grid;
      bh_ffr_build_resolution_radius_grid(p, &grid);
      const double theta_ref =
          BHFFRJetBroadenAnglesDeg[0] * M_PI / 180.0;
      const double omega_ref = 1.0 - cos(theta_ref);
      const double rref = grid.Radius[grid.Count - 1];

      double fplus = 0.0, fminus = 0.0;
      const int selected_active =
          bh_ffr_jet_select_live_broadening(res, &fplus, &fminus);

      /* If active targets are temporarily insufficient, keep a provisional
       * widest-angle geometry solely to define a threshold and queue receiver
       * wakeups. It is never allowed to fire until the live criterion passes. */
      const int selected =
          selected_active >= 0
              ? selected_active
              : BH_FFR_JET_BROADEN_ANGLE_COUNT - 1;

      const double theta =
          BHFFRJetBroadenAnglesDeg[selected] * M_PI / 180.0;
      const double radius =
          rref * cbrt(omega_ref / (1.0 - cos(theta)));

      if(!isfinite(radius) || !(radius > 0))
        terminate("BH_FFR: invalid selected live jet radius=%g angle=%g",
                  radius, BHFFRJetBroadenAnglesDeg[selected]);

      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      if(!isfinite(a) || !(a > 0))
        terminate("BH_FFR: invalid scale factor=%g in live jet setup", a);

      ev->AngleIndex = selected;
      ev->AngleDeg = BHFFRJetBroadenAnglesDeg[selected];
      ev->RadiusProper = radius;
      ev->RadiusCoordinate = radius / a;
      ev->CosCone = cos(theta);
      ev->Underresolved = selected_active < 0 ? 1 : 0;

      const double menc = res->BroadenEnclosedMass[selected];
      const double sigma2 = BHP[b].SigmaDM * BHP[b].SigmaDM;
      double vbind2 = sigma2;
      if(All.BHUseCentralBindingTerm)
        vbind2 +=
            2.0 * All.G * bh_ffr_central_mass_code(p) / radius;

      if(!isfinite(vbind2) || vbind2 < 0 ||
         !isfinite(menc) || menc < 0)
        terminate("BH_FFR: invalid live jet threshold environment "
                  "Menc=%g R=%g vbind2=%g for ID=%llu",
                  menc, radius, vbind2, (unsigned long long)P[p].ID);

      const double base_threshold = 0.5 * menc * vbind2;
      BHP[b].JetThresholdEnergy =
          All.BHJetBurstFactor * base_threshold;

      if(JetConflictRound == 0 && selected_active < 0)
        {
          printf("BH_FFR: jet live geometry ID=%llu task=%d "
                 "selectedAngle=none wakeAngle=%g Rwake=%g "
                 "underresolved=1 bufferRetained=1\n",
                 (unsigned long long)P[p].ID, ThisTask,
                 ev->AngleDeg, ev->RadiusProper);
          fflush(stdout);
        }

      if(JetConflictRound == 0 && BHP[b].SigmaDM > 0)
        {
          printf("BH_FFR: binding threshold jet ID=%llu task=%d "
                 "angle=%g Rjet=%g angleIndex=%d Menc=%g sigmaDM=%g "
                 "vbind2=%g central=%d EthJet=%g\n",
                 (unsigned long long)P[p].ID, ThisTask,
                 ev->AngleDeg, ev->RadiusProper, ev->AngleIndex,
                 menc, BHP[b].SigmaDM, vbind2,
                 All.BHUseCentralBindingTerm,
                 BHP[b].JetThresholdEnergy);
          fflush(stdout);
        }

      if(ev->Underresolved)
        continue;

      const double eth = BHP[b].JetThresholdEnergy;
      if(!(eth > 0) || !(BHP[b].JetEnergyBuffer >= eth))
        continue;

      ev->EnergyBefore = BHP[b].JetEnergyBuffer;
      ev->PacketEnergy = eth;

      for(int l = 0; l < 2; l++)
        {
          ev->LobeMass[l] = res->BroadenLobeMass[selected][l];
          ev->LobeCount[l] = res->BroadenLobeCount[selected][l];
        }

      const double pplus =
          res->BroadenLobeProjectedMomentum[selected][0];
      const double pminus =
          res->BroadenLobeProjectedMomentum[selected][1];

      if(!isfinite(ev->LobeMass[0]) || !isfinite(ev->LobeMass[1]) ||
         !(ev->LobeMass[0] > 0) || !(ev->LobeMass[1] > 0))
        terminate("BH_FFR: invalid selected live jet lobe masses ID=%llu "
                  "angle=%g Mplus=%g Mminus=%g",
                  (unsigned long long)P[p].ID, ev->AngleDeg,
                  (double)ev->LobeMass[0], (double)ev->LobeMass[1]);

      const double A =
          pplus / ev->LobeMass[0] - pminus / ev->LobeMass[1];
      const double B =
          1.0 / ev->LobeMass[0] + 1.0 / ev->LobeMass[1];

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
      printf("BH_FFR: jet coupling ID=%llu task=%d angle=%g Rjet=%g "
             "angleIndex=%d cosCone=%g Nplus=%lld Nminus=%lld\n",
             (unsigned long long)P[p].ID, ThisTask,
             (double)ev->AngleDeg, ev->RadiusProper, ev->AngleIndex,
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

  /* Constant-volume broadening regression:
   * R^3(1-cos theta) must be invariant across the survey angles. */
  const double rref = 64.0;
  const double theta_ref =
      BHFFRJetBroadenAnglesDeg[0] * M_PI / 180.0;
  const double vref =
      rref * rref * rref * (1.0 - cos(theta_ref));

  for(int qidx = 0; qidx < BH_FFR_JET_BROADEN_ANGLE_COUNT; qidx++)
    {
      const double theta =
          BHFFRJetBroadenAnglesDeg[qidx] * M_PI / 180.0;
      const double r =
          rref * cbrt((1.0 - cos(theta_ref)) / (1.0 - cos(theta)));
      const double v = r * r * r * (1.0 - cos(theta));
      if(!isfinite(r) || !(r > 0) ||
         fabs(v / vref - 1.0) > 3.0e-13)
        terminate("BH_FFR: jet broadening self-test failed q=%d theta=%g R=%g V/Vref=%g",
                  qidx, BHFFRJetBroadenAnglesDeg[qidx], r, v / vref);
    }
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

  for(int round = 0; round < BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP; round++)
    {
      JetConflictRound = round;
      bh_ffr_jet_comm_pass(BH_FFR_JET_STATS);
      bh_ffr_jet_report_broadening_survey();
      bh_ffr_jet_prepare_candidates();

      if(round == 0)
        bh_ffr_jet_comm_pass(BH_FFR_JET_WAKE);

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

  /* Mirror the wind drainage diagnostic for the jet reservoir. Residual
   * backlog is retained for later hydro updates once this step reaches the
   * bounded packet ceiling. */
  for(int n = 0; n < JetNTargets; n++)
    {
      const int p =
          bh_ffr_jet_particle_from_target(
              n, "bh_ffr_inject_jet_feedback/drain");
      const int b = P[p].BHDataIndex;
      const double eth = BHP[b].JetThresholdEnergy;
      const double ratio =
          eth > 0 ? BHP[b].JetEnergyBuffer / eth : 0.0;
      const int backlog = eth > 0 && BHP[b].JetEnergyBuffer >= eth;
      const int cap_hit =
          JetPacketCount[n] >= BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP;

      printf("BH_FFR: jet drain ID=%llu task=%d packets=%d maxPackets=%d "
             "buffer=%g Eth=%g backlogRatio=%g backlog=%d capHit=%d\n",
             (unsigned long long)P[p].ID, ThisTask,
             JetPacketCount[n], BH_FFR_FEEDBACK_MAX_PACKETS_PER_STEP,
             BHP[b].JetEnergyBuffer, eth, ratio, backlog, cap_hit);
      fflush(stdout);

    }

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
