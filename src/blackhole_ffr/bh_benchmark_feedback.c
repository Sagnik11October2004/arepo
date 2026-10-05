#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * IllustrisTNG thermal/kinetic feedback law for the convergence benchmark.
 *
 * Physics follows Weinberger et al. (2017):
 *
 *   chi = min[chi0 (M_BH/1e8 Msun)^beta, chi_max]
 *
 *   high state: Edot = eps_f,high eps_r Mdot_BH c^2
 *
 *   low state:  Edot = eps_f,kin Mdot_BH c^2
 *               eps_f,kin = min[rho/(f_thresh rho_SF), eps_f,kin,max]
 *
 *   E_inj,min = f_re (1/2) sigma_DM^2 M_enc
 *
 * and a kinetic event gives every target cell a kick in one common random
 * direction,
 *
 *   Delta p_j = m_j sqrt(2 DeltaE W_j/rho) n.
 *
 * The original TNG cosmological model adjusts its smoothing length to obtain a
 * prescribed neighbour count.  The convergence benchmark intentionally uses
 * the common fixed proper BHFeedbackRadius for every feedback prescription.
 * Inside that controlled aperture we retain TNG's cubic-spline weighting.
 *
 * Public AREPO does not provide the private TNG wake-up machinery.  As for the
 * validated FFR-MACER implementation, only synchronized active gas is modified.
 * The full aperture is read-only for M_enc and the density entering the TNG
 * efficiency.  If too little target mass is active, energy remains buffered.
 * Injection weights are renormalized over the active subset so no buffered
 * energy is silently lost.
 */

enum bh_tng_pass
{
  BH_TNG_PASS_STATS = 0,
  BH_TNG_PASS_MARK = 1,
  BH_TNG_PASS_CONFLICT = 2,
  BH_TNG_PASS_INJECT = 3
};

enum bh_tng_channel
{
  BH_TNG_CHANNEL_THERMAL = 0,
  BH_TNG_CHANNEL_KINETIC = 1
};

struct bh_tng_event
{
  int Candidate;
  int Fire;
  int Channel;
  MyDouble Energy;
  MyDouble Direction[3];
  MyDouble ActiveKernelDensity;
  MyDouble EnergyBefore;
};

typedef struct
{
  MyDouble Pos[3];
  MyFloat Radius;
  MyIDType BHID;
  int Candidate;
  int Fire;
  int Channel;
  MyDouble Energy;
  MyDouble Direction[3];
  MyDouble ActiveKernelDensity;
  int Firstnode;
} data_in;

typedef struct
{
  MyDouble EnclosedMass;
  MyDouble KernelDensityAll;
  MyDouble ActiveMass;
  MyDouble ActiveKernelDensity;
  long long ActiveCount;
  int Conflict;

  MyDouble InjectedThermalEnergy;
  MyDouble InjectedKineticQuadraticEnergy;
  MyDouble InjectedKineticLabEnergy;
  MyDouble InjectedMomentum[3];
} data_out;

static data_in *DataIn, *DataGet;
static data_out *DataResult, *DataOut;
static data_out *TNGResults;
static struct bh_tng_event *TNGEvents;
static MyIDType *TNGWinnerID;
static int TNGNTargets;
static int TNGPass;
static int TNGConflictRound;

static int bh_tng_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= TNGNTargets)
    terminate("BH_TNG: target=%d outside [0,%d) in %s", target, TNGNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_TNG: invalid active BH particle=%d in %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_TNG: invalid compact state for particle ID=%llu in %s",
              (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_TNG: particle ID=%llu not drifted to Ti_Current=%lld in %s",
              (unsigned long long)P[p].ID, (long long)All.Ti_Current, where);

  return p;
}

static int bh_tng_decode_hydro_timebin(int bin)
{
  if(bin < 0)
    bin = -bin - 1;
  return bin;
}

static int bh_tng_gas_is_active(int j)
{
  const int bin = bh_tng_decode_hydro_timebin(P[j].TimeBinHydro);
  return bin > 0 && bin < TIMEBINS && TimeBinSynchronized[bin];
}

static double bh_tng_coordinate_radius(void)
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0) || !isfinite(All.BHFeedbackRadius) ||
     !(All.BHFeedbackRadius > 0))
    terminate("BH_TNG: invalid feedback-radius conversion Rfb=%g a=%g",
              All.BHFeedbackRadius, a);

  const double r = All.BHFeedbackRadius / a;
  if(!isfinite(r) || !(r > 0))
    terminate("BH_TNG: invalid coordinate feedback radius=%g", r);
  return r;
}

/* Standard 3D cubic spline with support r<h, using the same normalization as
 * AREPO's KERNEL_COEFF_* constants.  Here r and h are proper code lengths. */
static double bh_tng_kernel(double r, double h)
{
  if(!isfinite(r) || r < 0 || !isfinite(h) || !(h > 0))
    terminate("BH_TNG: invalid kernel input r=%g h=%g", r, h);

  const double q = r / h;
  if(q >= 1.0)
    return 0.0;

  const double hinv = 1.0 / h;
  const double hinv3 = hinv * hinv * hinv;
  double w;

  if(q < 0.5)
    w = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (q - 1.0) * q * q);
  else
    {
      const double u = 1.0 - q;
      w = hinv3 * KERNEL_COEFF_5 * u * u * u;
    }

  if(!isfinite(w) || w < 0)
    terminate("BH_TNG: invalid kernel value W=%g q=%g", w, q);
  return w;
}

static double bh_tng_mode_threshold(double mass_msun, double chi0, double beta, double chimax)
{
  if(!isfinite(mass_msun) || mass_msun < 0 || !isfinite(chi0) || !(chi0 > 0) ||
     !isfinite(beta) || beta < 0 || !isfinite(chimax) || !(chimax > 0))
    terminate("BH_TNG: invalid mode-threshold inputs M=%g chi0=%g beta=%g chimax=%g",
              mass_msun, chi0, beta, chimax);

  double chi = chi0 * pow(mass_msun / 1.0e8, beta);
  if(chi > chimax)
    chi = chimax;
  if(!isfinite(chi) || chi < 0)
    terminate("BH_TNG: invalid mode threshold chi=%g", chi);
  return chi;
}

static double bh_tng_kinetic_efficiency(double nh, double density_factor,
                                        double sf_threshold_nh, double epsmax)
{
  if(!isfinite(nh) || nh < 0 || !isfinite(density_factor) || !(density_factor > 0) ||
     !isfinite(sf_threshold_nh) || !(sf_threshold_nh > 0) ||
     !isfinite(epsmax) || epsmax < 0)
    terminate("BH_TNG: invalid kinetic-efficiency inputs nH=%g f=%g nHSF=%g epsmax=%g",
              nh, density_factor, sf_threshold_nh, epsmax);

  const double eps = fmin(nh / (density_factor * sf_threshold_nh), epsmax);
  if(!isfinite(eps) || eps < 0)
    terminate("BH_TNG: invalid kinetic efficiency=%g", eps);
  return eps;
}

static uint64_t bh_tng_splitmix64(uint64_t x)
{
  x += UINT64_C(0x9e3779b97f4a7c15);
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

static double bh_tng_uniform01(uint64_t x)
{
  return (double)(bh_tng_splitmix64(x) >> 11) * (1.0 / 9007199254740992.0);
}

/* Reproducible pseudorandom realization of TNG's random per-event kick
 * direction.  This avoids dependence on MPI ordering or the global GSL RNG. */
static void bh_tng_random_direction(MyIDType id, integertime ti, int round, double n[3])
{
  uint64_t seed = (uint64_t)id;
  seed ^= bh_tng_splitmix64((uint64_t)ti + UINT64_C(0xd1b54a32d192ed03));
  seed ^= bh_tng_splitmix64((uint64_t)(round + 1) + UINT64_C(0x94d049bb133111eb));

  const double u = bh_tng_uniform01(seed);
  const double v = bh_tng_uniform01(seed ^ UINT64_C(0x632be59bd9b4e019));
  const double z = 2.0 * u - 1.0;
  const double phi = 2.0 * M_PI * v;
  const double rxy = sqrt(fmax(0.0, 1.0 - z * z));

  n[0] = rxy * cos(phi);
  n[1] = rxy * sin(phi);
  n[2] = z;

  const double norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if(!isfinite(norm) || fabs(norm - 1.0) > 2.0e-14)
    terminate("BH_TNG: invalid random direction norm=%g", norm);
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_tng_particle_from_target(target, "particle2in");
  const struct bh_tng_event *ev = &TNGEvents[target];

  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[p].Pos[k];
      in->Direction[k] = ev->Direction[k];
    }

  in->Radius = bh_tng_coordinate_radius();
  in->BHID = P[p].ID;
  in->Candidate = ev->Candidate;
  in->Fire = ev->Fire;
  in->Channel = ev->Channel;
  in->Energy = ev->Energy;
  in->ActiveKernelDensity = ev->ActiveKernelDensity;
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= TNGNTargets)
    terminate("BH_TNG: result target=%d outside [0,%d)", target, TNGNTargets);

  data_out *res = &TNGResults[target];

  if(mode == MODE_LOCAL_PARTICLES)
    *res = *out;
  else
    {
      res->EnclosedMass += out->EnclosedMass;
      res->KernelDensityAll += out->KernelDensityAll;
      res->ActiveMass += out->ActiveMass;
      res->ActiveKernelDensity += out->ActiveKernelDensity;
      res->ActiveCount += out->ActiveCount;
      if(out->Conflict)
        res->Conflict = 1;

      res->InjectedThermalEnergy += out->InjectedThermalEnergy;
      res->InjectedKineticQuadraticEnergy += out->InjectedKineticQuadraticEnergy;
      res->InjectedKineticLabEnergy += out->InjectedKineticLabEnergy;
      for(int k = 0; k < 3; k++)
        res->InjectedMomentum[k] += out->InjectedMomentum[k];
    }
}

#include "../utils/generic_comm_helpers2.h"

static int bh_tng_evaluate(int target, int mode, int threadid);

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
      if(target >= TNGNTargets)
        break;

      bh_tng_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_tng_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static unsigned long long bh_tng_feedback_priority(MyIDType id)
{
  uint64_t x = (uint64_t)id;
  x ^= (uint64_t)All.Ti_Current + UINT64_C(0x9e3779b97f4a7c15) +
       (uint64_t)(TNGConflictRound + 1) * UINT64_C(0xbf58476d1ce4e5b9);
  return (unsigned long long)bh_tng_splitmix64(x);
}

static int bh_tng_candidate_wins(MyIDType candidate, MyIDType incumbent)
{
  if(incumbent == 0)
    return 1;

  const unsigned long long pc = bh_tng_feedback_priority(candidate);
  const unsigned long long pi = bh_tng_feedback_priority(incumbent);
  return pc < pi || (pc == pi && candidate < incumbent);
}

static int bh_tng_evaluate(int target, int mode, int threadid)
{
  data_in local, *in;
  int numnodes, *firstnode;

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

  if(TNGPass != BH_TNG_PASS_STATS)
    {
      if((TNGPass == BH_TNG_PASS_INJECT && !in->Fire) ||
         ((TNGPass == BH_TNG_PASS_MARK || TNGPass == BH_TNG_PASS_CONFLICT) &&
          !in->Candidate))
        {
          if(mode == MODE_LOCAL_PARTICLES)
            out2particle(&out, target, MODE_LOCAL_PARTICLES);
          else
            DataResult[target] = out;
          return 0;
        }
    }

  const int nfound =
      ngb_treefind_variable_threads(in->Pos, in->Radius, target, mode, threadid,
                                    numnodes, firstnode);

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_TNG: invalid scale factor=%g", a);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_TNG: gas index=%d outside NumGas=%d", j, NumGas);
      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      double xtmp, ytmp, ztmp;
      const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
      const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
      const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
      const double r2 = dx * dx + dy * dy + dz * dz;
      if(r2 > in->Radius * in->Radius)
        continue;

      const double rproper = sqrt(r2) * a;
      const double w = bh_tng_kernel(rproper, All.BHFeedbackRadius);
      if(!(w > 0))
        continue;

      const int active = bh_tng_gas_is_active(j);

      if(TNGPass == BH_TNG_PASS_STATS)
        {
          out.EnclosedMass += P[j].Mass;
          out.KernelDensityAll += P[j].Mass * w;
          if(active)
            {
              out.ActiveMass += P[j].Mass;
              out.ActiveKernelDensity += P[j].Mass * w;
              out.ActiveCount++;
            }
          continue;
        }

      if(!active)
        continue;

      if(TNGPass == BH_TNG_PASS_MARK)
        {
          if(bh_tng_candidate_wins(in->BHID, TNGWinnerID[j]))
            TNGWinnerID[j] = in->BHID;
          continue;
        }

      if(TNGPass == BH_TNG_PASS_CONFLICT)
        {
          if(TNGWinnerID[j] != in->BHID)
            out.Conflict = 1;
          continue;
        }

      if(TNGPass != BH_TNG_PASS_INJECT)
        terminate("BH_TNG: unknown pass=%d", TNGPass);

      if(!(in->ActiveKernelDensity > 0) || !(in->Energy > 0))
        terminate("BH_TNG: invalid injection normalization rhoActive=%g E=%g",
                  in->ActiveKernelDensity, in->Energy);

      if(in->Channel == BH_TNG_CHANNEL_THERMAL)
        {
          const double weight = P[j].Mass * w / in->ActiveKernelDensity;
          const double denergy = in->Energy * weight;
          if(!isfinite(weight) || weight < 0 || !isfinite(denergy) || denergy < 0)
            terminate("BH_TNG: invalid thermal weight=%g dE=%g for gas ID=%llu",
                      weight, denergy, (unsigned long long)P[j].ID);

          SphP[j].Energy += a * a * denergy;
          out.InjectedThermalEnergy += denergy;
        }
      else if(in->Channel == BH_TNG_CHANNEL_KINETIC)
        {
          const double dvkick =
              sqrt(2.0 * in->Energy * w / in->ActiveKernelDensity);
          if(!isfinite(dvkick) || dvkick < 0)
            terminate("BH_TNG: invalid kinetic kick speed=%g for gas ID=%llu",
                      dvkick, (unsigned long long)P[j].ID);

          double dp_phys[3], v_phys[3];
          double dp2 = 0.0, cross = 0.0;
          for(int k = 0; k < 3; k++)
            {
              v_phys[k] = (SphP[j].Momentum[k] / P[j].Mass) / a;
              dp_phys[k] = P[j].Mass * dvkick * in->Direction[k];
              dp2 += dp_phys[k] * dp_phys[k];
              cross += v_phys[k] * dp_phys[k];
            }

          const double quadratic = 0.5 * dp2 / P[j].Mass;
          const double dkin_lab = cross + quadratic;
          if(!isfinite(quadratic) || quadratic < 0 || !isfinite(dkin_lab))
            terminate("BH_TNG: invalid kinetic increments quadratic=%g lab=%g for gas ID=%llu",
                      quadratic, dkin_lab, (unsigned long long)P[j].ID);

          for(int k = 0; k < 3; k++)
            {
              SphP[j].Momentum[k] += a * dp_phys[k];
              out.InjectedMomentum[k] += dp_phys[k];
            }
          SphP[j].Energy += a * a * dkin_lab;

          for(int k = 0; k < 3; k++)
            P[j].Vel[k] = SphP[j].Momentum[k] / P[j].Mass;

          out.InjectedKineticQuadraticEnergy += quadratic;
          out.InjectedKineticLabEnergy += dkin_lab;
        }
      else
        terminate("BH_TNG: unknown feedback channel=%d", in->Channel);
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static void bh_tng_comm_pass(int pass)
{
  TNGPass = pass;
  memset(TNGResults, 0, (TNGNTargets > 0 ? TNGNTargets : 1) * sizeof(*TNGResults));

  if(All.TotNumGas > 0)
    {
      generic_set_MaxNexport();
      generic_comm_pattern(TNGNTargets, kernel_local, kernel_imported);
    }
}

static void bh_tng_allocate_work(void)
{
  TNGNTargets = NumActiveBHFFR;
  TNGEvents = (struct bh_tng_event *)mymalloc(
      "BHTNGEvents", (TNGNTargets > 0 ? TNGNTargets : 1) * sizeof(*TNGEvents));
  TNGResults = (data_out *)mymalloc(
      "BHTNGResults", (TNGNTargets > 0 ? TNGNTargets : 1) * sizeof(*TNGResults));
  TNGWinnerID = (MyIDType *)mymalloc(
      "BHTNGWinner", (NumGas > 0 ? NumGas : 1) * sizeof(*TNGWinnerID));

  memset(TNGEvents, 0, (TNGNTargets > 0 ? TNGNTargets : 1) * sizeof(*TNGEvents));
  memset(TNGResults, 0, (TNGNTargets > 0 ? TNGNTargets : 1) * sizeof(*TNGResults));
  memset(TNGWinnerID, 0, (NumGas > 0 ? NumGas : 1) * sizeof(*TNGWinnerID));
}

static void bh_tng_free_work(void)
{
  myfree(TNGWinnerID);
  myfree(TNGResults);
  myfree(TNGEvents);
  TNGWinnerID = NULL;
  TNGResults = NULL;
  TNGEvents = NULL;
  TNGNTargets = 0;
}

static double bh_tng_density_to_nh(double density_code)
{
  if(!isfinite(density_code) || density_code < 0)
    terminate("BH_TNG: invalid density=%g", density_code);

  const double rho_cgs =
      density_code * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  const double nh = HYDROGEN_MASSFRAC * rho_cgs / PROTONMASS;
  if(!isfinite(nh) || nh < 0)
    terminate("BH_TNG: invalid hydrogen number density=%g", nh);
  return nh;
}

void bh_benchmark_tng_feedback_accumulate(void)
{
  if(All.BHBenchmarkFeedbackModel != BH_BENCHMARK_FEEDBACK_TNG)
    return;

  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

  bh_tng_allocate_work();
  bh_tng_comm_pass(BH_TNG_PASS_STATS);

  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  if(!isfinite(c_internal) || !(c_internal > 0))
    terminate("BH_TNG: invalid internal light speed=%g", c_internal);

  for(int n = 0; n < TNGNTargets; n++)
    {
      const int p = bh_tng_particle_from_target(n, "feedback_accumulate");
      const int b = P[p].BHDataIndex;
      const data_out *res = &TNGResults[n];

      const double bh_mass_msun =
          BHP[b].BHMass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
      const double chi =
          bh_tng_mode_threshold(bh_mass_msun, All.BHBenchmarkTNGChi0,
                                All.BHBenchmarkTNGChiBeta,
                                All.BHBenchmarkTNGChiMax);
      const double fedd_raw =
          BHP[b].MdotEddington > 0
              ? BHP[b].BenchmarkMdotRaw / BHP[b].MdotEddington
              : 0.0;

      if(!isfinite(fedd_raw) || fedd_raw < 0)
        terminate("BH_TNG: invalid raw Eddington ratio=%g for ID=%llu",
                  fedd_raw, (unsigned long long)P[p].ID);

      const double nh = bh_tng_density_to_nh(res->KernelDensityAll);
      const double epskin =
          bh_tng_kinetic_efficiency(nh,
                                    All.BHBenchmarkTNGKineticDensityFactor,
                                    All.BHBenchmarkTNGSFThresholdNH,
                                    All.BHBenchmarkTNGKineticMaxEfficiency);

      const int mode =
          fedd_raw >= chi ? BH_BENCHMARK_TNG_MODE_THERMAL
                          : BH_BENCHMARK_TNG_MODE_KINETIC;

      double power;
      if(mode == BH_BENCHMARK_TNG_MODE_THERMAL)
        power = All.BHBenchmarkTNGThermalCoupling *
                All.BHBenchmarkRadiativeEfficiency *
                BHP[b].BenchmarkMdotOperational * c_internal * c_internal;
      else
        power = epskin * BHP[b].BenchmarkMdotOperational * c_internal * c_internal;

      if(!isfinite(power) || power < 0)
        terminate("BH_TNG: invalid feedback power=%g for ID=%llu",
                  power, (unsigned long long)P[p].ID);

      double dt_code = 0.0;
      if(BHP[b].LastProcessedTi != All.Ti_Current)
        dt_code = bh_ffr_get_elapsed_time_code_time(p);
      if(!isfinite(dt_code) || dt_code < 0)
        terminate("BH_TNG: invalid elapsed code time=%g for ID=%llu",
                  dt_code, (unsigned long long)P[p].ID);

      const double denergy = power * dt_code;
      if(!isfinite(denergy) || denergy < 0)
        terminate("BH_TNG: invalid accumulated energy=%g for ID=%llu",
                  denergy, (unsigned long long)P[p].ID);

      if(mode == BH_BENCHMARK_TNG_MODE_THERMAL)
        BHP[b].TNGThermalEnergyBuffer += denergy;
      else
        BHP[b].TNGKineticEnergyBuffer += denergy;

      const double eth =
          0.5 * All.BHBenchmarkTNGKineticBurstFactor * res->EnclosedMass *
          BHP[b].SigmaDM * BHP[b].SigmaDM;

      if(!isfinite(eth) || eth < 0)
        terminate("BH_TNG: invalid kinetic threshold=%g for ID=%llu",
                  eth, (unsigned long long)P[p].ID);

      BHP[b].TNGFeedbackPower = power;
      BHP[b].TNGKineticThresholdEnergy = eth;
      BHP[b].TNGEddingtonRatio = fedd_raw;
      BHP[b].TNGModeThreshold = chi;
      BHP[b].TNGFeedbackMode = mode;

      printf("BH_TNG: source ID=%llu task=%d mode=%s MbhMsun=%g "
             "fEddRaw=%g chi=%g mdotRaw=%g mdotOperational=%g mdotRealized=%g "
             "nH=%g epsKin=%g power=%g dE=%g Eth=%g EthermBuf=%g EkinBuf=%g "
             "Menc=%g sigmaDM=%g\n",
             (unsigned long long)P[p].ID, ThisTask,
             mode == BH_BENCHMARK_TNG_MODE_THERMAL ? "thermal" : "kinetic",
             bh_mass_msun, fedd_raw, chi, BHP[b].BenchmarkMdotRaw,
             BHP[b].BenchmarkMdotOperational, BHP[b].MdotHorizon, nh, epskin,
             power, denergy, eth,
             BHP[b].TNGThermalEnergyBuffer, BHP[b].TNGKineticEnergyBuffer,
             res->EnclosedMass, BHP[b].SigmaDM);
      fflush(stdout);
    }

  bh_tng_free_work();
}

static int bh_tng_global_candidate_count(void)
{
  int local = 0;
  for(int n = 0; n < TNGNTargets; n++)
    if(TNGEvents[n].Candidate)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static int bh_tng_global_fire_count(void)
{
  int local = 0;
  for(int n = 0; n < TNGNTargets; n++)
    if(TNGEvents[n].Fire)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static void bh_tng_prepare_events(int channel, int round)
{
  for(int n = 0; n < TNGNTargets; n++)
    {
      struct bh_tng_event *ev = &TNGEvents[n];
      memset(ev, 0, sizeof(*ev));

      const int p = bh_tng_particle_from_target(n, "prepare_events");
      const int b = P[p].BHDataIndex;
      const data_out *res = &TNGResults[n];

      ev->Channel = channel;
      ev->ActiveKernelDensity = res->ActiveKernelDensity;

      if(res->ActiveCount <= 0 || !(res->ActiveKernelDensity > 0) ||
         res->ActiveMass < All.BHMinActiveTargetMassFrac * res->EnclosedMass)
        continue;

      if(channel == BH_TNG_CHANNEL_THERMAL)
        {
          if(!(BHP[b].TNGThermalEnergyBuffer > 0))
            continue;

          ev->EnergyBefore = BHP[b].TNGThermalEnergyBuffer;
          ev->Energy = ev->EnergyBefore;
          ev->Candidate = 1;
        }
      else
        {
          const double eth =
              0.5 * All.BHBenchmarkTNGKineticBurstFactor * res->EnclosedMass *
              BHP[b].SigmaDM * BHP[b].SigmaDM;
          BHP[b].TNGKineticThresholdEnergy = eth;

          if(!(eth > 0) || BHP[b].TNGKineticEnergyBuffer < eth)
            continue;

          ev->EnergyBefore = BHP[b].TNGKineticEnergyBuffer;
          ev->Energy = ev->EnergyBefore;
          bh_tng_random_direction(P[p].ID, All.Ti_Current, round, ev->Direction);
          ev->Candidate = 1;
        }
    }
}

static void bh_tng_resolve_conflicts(void)
{
  const int candidates = bh_tng_global_candidate_count();
  if(candidates <= 0)
    return;

  if(candidates == 1)
    {
      for(int n = 0; n < TNGNTargets; n++)
        TNGEvents[n].Fire = TNGEvents[n].Candidate;
      return;
    }

  memset(TNGWinnerID, 0, (NumGas > 0 ? NumGas : 1) * sizeof(*TNGWinnerID));
  bh_tng_comm_pass(BH_TNG_PASS_MARK);
  bh_tng_comm_pass(BH_TNG_PASS_CONFLICT);

  for(int n = 0; n < TNGNTargets; n++)
    TNGEvents[n].Fire = TNGEvents[n].Candidate && !TNGResults[n].Conflict;
}

static void bh_tng_commit_events(int channel)
{
  for(int n = 0; n < TNGNTargets; n++)
    {
      const struct bh_tng_event *ev = &TNGEvents[n];
      if(!ev->Fire)
        continue;

      const int p = bh_tng_particle_from_target(n, "commit_events");
      const int b = P[p].BHDataIndex;
      const data_out *res = &TNGResults[n];

      if(channel == BH_TNG_CHANNEL_THERMAL)
        {
          const double tol = 3.0e-8 * fmax(fabs(ev->Energy), 1.0e-30);
          if(fabs(res->InjectedThermalEnergy - ev->Energy) > tol)
            terminate("BH_TNG: thermal energy mismatch ID=%llu injected=%g requested=%g",
                      (unsigned long long)P[p].ID,
                      res->InjectedThermalEnergy, ev->Energy);

          BHP[b].TNGThermalEnergyBuffer = ev->EnergyBefore - ev->Energy;
          if(BHP[b].TNGThermalEnergyBuffer < 0 &&
             BHP[b].TNGThermalEnergyBuffer >
                 -3.0e-12 * fmax(fabs(ev->EnergyBefore), 1.0e-30))
            BHP[b].TNGThermalEnergyBuffer = 0.0;
          if(BHP[b].TNGThermalEnergyBuffer < 0)
            terminate("BH_TNG: negative thermal buffer after event for ID=%llu",
                      (unsigned long long)P[p].ID);

          printf("BH_TNG: thermal event ID=%llu task=%d E=%g dErel=%g "
                 "Nactive=%lld Mactive=%g\n",
                 (unsigned long long)P[p].ID, ThisTask, ev->Energy,
                 fabs(res->InjectedThermalEnergy - ev->Energy) /
                     fmax(ev->Energy, 1.0e-30),
                 TNGResults[n].ActiveCount, TNGResults[n].ActiveMass);
          fflush(stdout);
        }
      else
        {
          /* Equation (10) guarantees that the sum of the quadratic kinetic
           * terms is DeltaE.  The lab-frame total-energy change contains the
           * pre-existing-velocity cross term and therefore is not generally
           * DeltaE for a single TNG event. */
          const double tol = 3.0e-8 * fmax(fabs(ev->Energy), 1.0e-30);
          if(fabs(res->InjectedKineticQuadraticEnergy - ev->Energy) > tol)
            terminate("BH_TNG: kinetic quadratic-energy mismatch ID=%llu injected=%g requested=%g",
                      (unsigned long long)P[p].ID,
                      res->InjectedKineticQuadraticEnergy, ev->Energy);

          BHP[b].TNGKineticEnergyBuffer = ev->EnergyBefore - ev->Energy;
          if(BHP[b].TNGKineticEnergyBuffer < 0 &&
             BHP[b].TNGKineticEnergyBuffer >
                 -3.0e-12 * fmax(fabs(ev->EnergyBefore), 1.0e-30))
            BHP[b].TNGKineticEnergyBuffer = 0.0;
          if(BHP[b].TNGKineticEnergyBuffer < 0)
            terminate("BH_TNG: negative kinetic buffer after event for ID=%llu",
                      (unsigned long long)P[p].ID);

          double pnorm2 = 0.0;
          for(int k = 0; k < 3; k++)
            pnorm2 += res->InjectedMomentum[k] * res->InjectedMomentum[k];

          printf("BH_TNG: kinetic event ID=%llu task=%d E=%g dEquadRel=%g "
                 "dElab=%g px=%g py=%g pz=%g pnorm=%g "
                 "dir=(%g,%g,%g)\n",
                 (unsigned long long)P[p].ID, ThisTask, ev->Energy,
                 fabs(res->InjectedKineticQuadraticEnergy - ev->Energy) /
                     fmax(ev->Energy, 1.0e-30),
                 res->InjectedKineticLabEnergy,
                 res->InjectedMomentum[0], res->InjectedMomentum[1],
                 res->InjectedMomentum[2], sqrt(pnorm2),
                 ev->Direction[0], ev->Direction[1], ev->Direction[2]);
          fflush(stdout);
        }
    }
}

static int bh_tng_process_channel(int channel, int round)
{
  bh_tng_comm_pass(BH_TNG_PASS_STATS);
  bh_tng_prepare_events(channel, round);

  if(bh_tng_global_candidate_count() <= 0)
    return 0;

  TNGConflictRound = round;
  bh_tng_resolve_conflicts();
  const int fires = bh_tng_global_fire_count();
  if(fires <= 0)
    return 0;

  bh_tng_comm_pass(BH_TNG_PASS_INJECT);
  bh_tng_commit_events(channel);
  return fires;
}

void bh_benchmark_tng_feedback_inject(void)
{
  if(All.BHBenchmarkFeedbackModel != BH_BENCHMARK_FEEDBACK_TNG)
    return;

  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

  bh_tng_allocate_work();

  int fires = 0;
  fires += bh_tng_process_channel(BH_TNG_CHANNEL_THERMAL, 0);
  fires += bh_tng_process_channel(BH_TNG_CHANNEL_KINETIC, 1);

  if(fires > 0)
    update_primitive_variables();

  bh_ffr_validate_state("post-TNG-feedback");
  bh_tng_free_work();
}

void bh_benchmark_tng_feedback_self_test(void)
{
  const double chi_pivot = bh_tng_mode_threshold(1.0e8, 0.002, 2.0, 0.1);
  if(fabs(chi_pivot - 0.002) > 2.0e-14)
    terminate("BH_TNG: self-test chi pivot got=%g expected=0.002", chi_pivot);

  const double chi_cap = bh_tng_mode_threshold(1.0e10, 0.002, 2.0, 0.1);
  if(fabs(chi_cap - 0.1) > 2.0e-14)
    terminate("BH_TNG: self-test chi cap got=%g expected=0.1", chi_cap);

  const double eps_low = bh_tng_kinetic_efficiency(0.0005, 0.05, 0.1, 0.2);
  if(fabs(eps_low - 0.1) > 2.0e-14)
    terminate("BH_TNG: self-test low-density efficiency got=%g expected=0.1", eps_low);

  const double eps_cap = bh_tng_kinetic_efficiency(1.0, 0.05, 0.1, 0.2);
  if(fabs(eps_cap - 0.2) > 2.0e-14)
    terminate("BH_TNG: self-test efficiency cap got=%g expected=0.2", eps_cap);

  double dir[3];
  bh_tng_random_direction((MyIDType)1234567, (integertime)42, 0, dir);
  const double norm =
      sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  if(fabs(norm - 1.0) > 2.0e-14)
    terminate("BH_TNG: self-test random direction norm=%g", norm);

  /* Synthetic two-cell check of equation (10): rho=sum mW implies that the
   * sum 1/2 dp^2/m is exactly the requested event energy. */
  const double m1 = 2.0, m2 = 3.0, w1 = 0.7, w2 = 0.2, e = 5.0;
  const double rho = m1 * w1 + m2 * w2;
  const double dv1 = sqrt(2.0 * e * w1 / rho);
  const double dv2 = sqrt(2.0 * e * w2 / rho);
  const double got = 0.5 * m1 * dv1 * dv1 + 0.5 * m2 * dv2 * dv2;
  if(fabs(got - e) > 2.0e-13 * e)
    terminate("BH_TNG: self-test kinetic normalization got=%g expected=%g",
              got, e);
}
