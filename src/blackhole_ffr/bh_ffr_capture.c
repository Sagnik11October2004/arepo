#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/* Three-pass benchmark capture:
 * (1) a common fixed-aperture environment, (2) overlap-safe sink coefficients,
 * (3) conservative per-BH shares.  The validated FFR path is recovered when
 * BHBenchmarkAccretionModel=FFR and the target is the reservoir. */

enum bh_ffr_capture_pass
{
  BH_FFR_CAPTURE_PASS_ENVIRONMENT = 0,
  BH_FFR_CAPTURE_PASS_LAMBDA = 1,
  BH_FFR_CAPTURE_PASS_SHARES = 2
};

struct bh_ffr_capture_result
{
  MyDouble CapturedMass;
  MyDouble CapturedRate;
  MyDouble CapturedMomentum[3];
  MyDouble CapturedAbsMomentum[3];
  MyDouble CapturedKineticEnergy;
  MyDouble CoherenceIncrement[3];

  /* Transient common-aperture environment used by every benchmark model. */
  MyDouble EnvGasMass;
  MyDouble EnvVolume;
  MyDouble EnvSoundVolumeWeighted;
  MyDouble EnvVelocityVolumeWeighted[3];
  MyDouble EnvAngularMomentum[3];
  MyDouble EnvFFRRawRate;

  /* Algebraic model result before the conservative cell sink is applied. */
  MyDouble ModelRawRate;
  MyDouble ModelOperationalRate;
  MyDouble ModelEddingtonRate;
  MyDouble ModelBoostFactor;
  MyDouble ModelAMLimiter;
  MyDouble LambdaScale;
  MyDouble UniformLambda;

  int MinHydroTimeBin;
};

static double *LambdaSink;
static double *SinkMass;
static struct bh_ffr_capture_result *CaptureResults;
static int CaptureNTargets;
static int CapturePass;
static int CaptureSelfTestDone;

typedef struct
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat AccretionRadius;
  MyDouble CentralMass;
  MyDouble SchwarzschildRadius;
  MyDouble LambdaScale;
  MyDouble UniformLambda;
  int CaptureEnabled;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble CapturedMass;
  MyDouble CapturedRate;
  MyDouble CapturedMomentum[3];
  MyDouble CapturedAbsMomentum[3];
  MyDouble CapturedKineticEnergy;
  MyDouble CoherenceIncrement[3];
  MyDouble EnvGasMass;
  MyDouble EnvVolume;
  MyDouble EnvSoundVolumeWeighted;
  MyDouble EnvVelocityVolumeWeighted[3];
  MyDouble EnvAngularMomentum[3];
  MyDouble EnvFFRRawRate;
  int MinHydroTimeBin;
} data_out;

static data_out *DataResult, *DataOut;

static int bh_ffr_decode_hydro_timebin(int bin)
{
  if(bin < 0)
    bin = -bin - 1;
  return bin;
}

static int bh_ffr_gas_is_active(int j)
{
  const int bin = bh_ffr_decode_hydro_timebin(P[j].TimeBinHydro);
  return bin > 0 && bin < TIMEBINS && TimeBinSynchronized[bin];
}

static double bh_ffr_proper_radius_to_coordinate_radius(double proper_radius)
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;

  if(!isfinite(a) || a <= 0 || !isfinite(proper_radius) || proper_radius <= 0)
    terminate("BH_FFR: invalid proper-radius conversion R=%g a=%g", proper_radius, a);

  return proper_radius / a;
}

static int bh_ffr_capture_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= CaptureNTargets)
    terminate("BH_FFR: capture target %d outside [0,%d) in %s", target, CaptureNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active capture particle %d in %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for capture particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: capture particle ID=%llu is not drifted to Ti_Current=%lld in %s", (unsigned long long)P[p].ID,
              (long long)All.Ti_Current, where);

  return p;
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_capture_particle_from_target(target, "particle2in");
  const int b = P[p].BHDataIndex;

  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[p].Pos[k];
      in->Vel[k] = P[p].Vel[k];
    }

  in->AccretionRadius = bh_ffr_proper_radius_to_coordinate_radius(All.BHAccretionRadius);
  /* Expose the PDF source-fidelity M_cen=M_BH comparison while keeping
   * the full unresolved dynamical point mass as the recommended default. */
  in->CentralMass = bh_ffr_central_mass_code(p);

  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  in->SchwarzschildRadius =
      (BHP[b].BHMass > 0 && c_internal > 0) ? 2.0 * All.G * BHP[b].BHMass / (c_internal * c_internal) : 0.0;

  in->CaptureEnabled = (BHP[b].LastProcessedTi < All.Ti_Current);
  in->LambdaScale = 1.0;
  in->UniformLambda = 0.0;
  if(CapturePass != BH_FFR_CAPTURE_PASS_ENVIRONMENT)
    {
      if(target < 0 || target >= CaptureNTargets)
        terminate("BH_BENCHMARK: invalid capture target=%d while loading sink coefficients", target);
      in->LambdaScale = CaptureResults[target].LambdaScale;
      in->UniformLambda = CaptureResults[target].UniformLambda;
    }
  in->Firstnode = firstnode;

  if(!isfinite(in->CentralMass) || in->CentralMass < 0 || !isfinite(in->SchwarzschildRadius) ||
     in->SchwarzschildRadius < 0)
    terminate("BH_FFR: invalid capture mass/radius for particle ID=%llu", (unsigned long long)P[p].ID);
}

static void out2particle(data_out *out, int target, int mode)
{
  if(CapturePass == BH_FFR_CAPTURE_PASS_LAMBDA)
    return;

  if(target < 0 || target >= CaptureNTargets)
    terminate("BH_FFR: capture result target %d outside [0,%d)", target, CaptureNTargets);

  struct bh_ffr_capture_result *res = &CaptureResults[target];

  if(CapturePass == BH_FFR_CAPTURE_PASS_ENVIRONMENT)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          res->EnvGasMass = out->EnvGasMass;
          res->EnvVolume = out->EnvVolume;
          res->EnvSoundVolumeWeighted = out->EnvSoundVolumeWeighted;
          res->EnvFFRRawRate = out->EnvFFRRawRate;
          for(int k = 0; k < 3; k++)
            {
              res->EnvVelocityVolumeWeighted[k] = out->EnvVelocityVolumeWeighted[k];
              res->EnvAngularMomentum[k] = out->EnvAngularMomentum[k];
            }
        }
      else
        {
          res->EnvGasMass += out->EnvGasMass;
          res->EnvVolume += out->EnvVolume;
          res->EnvSoundVolumeWeighted += out->EnvSoundVolumeWeighted;
          res->EnvFFRRawRate += out->EnvFFRRawRate;
          for(int k = 0; k < 3; k++)
            {
              res->EnvVelocityVolumeWeighted[k] += out->EnvVelocityVolumeWeighted[k];
              res->EnvAngularMomentum[k] += out->EnvAngularMomentum[k];
            }
        }
      return;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    {
      res->CapturedMass = out->CapturedMass;
      res->CapturedRate = out->CapturedRate;
      res->CapturedKineticEnergy = out->CapturedKineticEnergy;
      for(int k = 0; k < 3; k++)
        {
          res->CapturedMomentum[k] = out->CapturedMomentum[k];
          res->CapturedAbsMomentum[k] = out->CapturedAbsMomentum[k];
          res->CoherenceIncrement[k] = out->CoherenceIncrement[k];
        }
      res->MinHydroTimeBin = out->MinHydroTimeBin;
    }
  else
    {
      res->CapturedMass += out->CapturedMass;
      res->CapturedRate += out->CapturedRate;
      res->CapturedKineticEnergy += out->CapturedKineticEnergy;
      for(int k = 0; k < 3; k++)
        {
          res->CapturedMomentum[k] += out->CapturedMomentum[k];
          res->CapturedAbsMomentum[k] += out->CapturedAbsMomentum[k];
          res->CoherenceIncrement[k] += out->CoherenceIncrement[k];
        }
      if(out->MinHydroTimeBin < res->MinHydroTimeBin)
        res->MinHydroTimeBin = out->MinHydroTimeBin;
    }
}

#include "../utils/generic_comm_helpers2.h"

static double bh_ffr_active_gas_timestep_code_time(int j);
static int bh_ffr_capture_evaluate(int target, int mode, int threadid);

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
      if(target >= CaptureNTargets)
        break;

      bh_ffr_capture_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_capture_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static double bh_ffr_capture_lambda_core(const data_in *bh, double coordinate_distance,
                                             double freefall_a, double freefall_alpha)
{
  if(!isfinite(freefall_a) || freefall_a < 0 || !isfinite(freefall_alpha))
    terminate("BH_FFR: invalid free-fall coefficient inputs A=%g alpha=%g", freefall_a, freefall_alpha);

  if(!bh->CaptureEnabled || freefall_a == 0 || !(bh->CentralMass > 0))
    return 0.0;

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  double d = coordinate_distance * a;

  const double d_floor = 1.0e-12 * All.BHAccretionRadius;
  if(d < d_floor)
    d = d_floor;

  const double tff = sqrt(d * d * d / (All.G * bh->CentralMass));
  if(!isfinite(tff) || !(tff > 0))
    terminate("BH_FFR: invalid free-fall time d=%g Mcen=%g G=%g", d, bh->CentralMass, All.G);

  double eta = freefall_a;
  if(freefall_alpha != 0)
    {
      if(!(bh->SchwarzschildRadius > 0))
        return 0.0;

      const double ratio = d / bh->SchwarzschildRadius;
      if(!(ratio > 0) || !isfinite(ratio))
        terminate("BH_FFR: invalid d/Rs=%g in free-fall sink", ratio);

      eta *= pow(ratio, freefall_alpha);
    }

  const double lambda = eta / tff;
  if(!isfinite(lambda) || lambda < 0)
    terminate("BH_FFR: invalid free-fall lambda=%g eta=%g tff=%g", lambda, eta, tff);

  return lambda;
}

static double bh_ffr_capture_lambda(const data_in *bh, double coordinate_distance)
{
  if(!bh->CaptureEnabled)
    return 0.0;

  if(All.BHBenchmarkAccretionModel == BH_BENCHMARK_ACC_FFR)
    return bh->LambdaScale *
           bh_ffr_capture_lambda_core(bh, coordinate_distance, All.BHFreeFallA, All.BHFreeFallAlpha);

  if(!isfinite(bh->UniformLambda) || bh->UniformLambda < 0)
    terminate("BH_BENCHMARK: invalid uniform sink lambda=%g", bh->UniformLambda);
  return bh->UniformLambda;
}

static int bh_ffr_capture_evaluate(int target, int mode, int threadid)
{
  data_in local, *bh;
  int numnodes, *firstnode;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      bh = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      bh = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  data_out out;
  memset(&out, 0, sizeof(out));
  out.MinHydroTimeBin = TIMEBINS;

  const int nfound =
      ngb_treefind_variable_threads(bh->Pos, bh->AccretionRadius, target, mode, threadid, numnodes, firstnode);

  if(nfound < 0)
    terminate("BH_FFR: neighbour search failed during capture pass %d", CapturePass);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_FFR: capture gas index %d outside NumGas=%d", j, NumGas);

      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      const int bin = bh_ffr_decode_hydro_timebin(P[j].TimeBinHydro);
      const double r = sqrt(Thread[threadid].R2list[n]);

      if(CapturePass == BH_FFR_CAPTURE_PASS_ENVIRONMENT)
        {
          if(!(SphP[j].Volume > 0) || !isfinite(SphP[j].Volume))
            terminate("BH_BENCHMARK: invalid gas volume=%g for ID=%llu", (double)SphP[j].Volume,
                      (unsigned long long)P[j].ID);

          const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
          const double cs = get_sound_speed(j);
          if(!isfinite(cs) || cs < 0)
            terminate("BH_BENCHMARK: invalid sound speed=%g for gas ID=%llu", cs, (unsigned long long)P[j].ID);

          out.EnvGasMass += P[j].Mass;
          out.EnvVolume += SphP[j].Volume;
          out.EnvSoundVolumeWeighted += SphP[j].Volume * cs;

          /* AREPO's NEAREST_[XYZ] periodic-wrapping macros use these
           * scratch variables internally, so they must exist in scope. */
          double xtmp, ytmp, ztmp;
          double dr[3] = {NEAREST_X(P[j].Pos[0] - bh->Pos[0]) * a,
                          NEAREST_Y(P[j].Pos[1] - bh->Pos[1]) * a,
                          NEAREST_Z(P[j].Pos[2] - bh->Pos[2]) * a};
          double dv[3];
          for(int k = 0; k < 3; k++)
            {
              const double vphys = P[j].Vel[k] / a;
              const double bhvphys = bh->Vel[k] / a;
              out.EnvVelocityVolumeWeighted[k] += SphP[j].Volume * vphys;
              dv[k] = vphys - bhvphys;
            }

          const double ell[3] = {dr[1] * dv[2] - dr[2] * dv[1],
                                 dr[2] * dv[0] - dr[0] * dv[2],
                                 dr[0] * dv[1] - dr[1] * dv[0]};
          for(int k = 0; k < 3; k++)
            out.EnvAngularMomentum[k] += P[j].Mass * ell[k];

          /* The benchmark raw rate is an instantaneous estimator and
           * should be defined even at the initial synchronization point, where
           * no elapsed transaction exists yet.  Keep CaptureEnabled gating for
           * the actual sink passes, but remove it from this diagnostic-only
           * environment estimate so FFR is comparable to the Bondi-family
           * raw rates at identical states. */
          data_in rate_bh = *bh;
          rate_bh.CaptureEnabled = 1;
          const double lambda_ffr =
              bh_ffr_capture_lambda_core(&rate_bh, r, All.BHFreeFallA, All.BHFreeFallAlpha);
          if(lambda_ffr > 0)
            out.EnvFFRRawRate += P[j].Mass * lambda_ffr;

          continue;
        }

      if(CapturePass == BH_FFR_CAPTURE_PASS_SHARES && bin > 0 && bin < out.MinHydroTimeBin)
        out.MinHydroTimeBin = bin;

      if(!bh_ffr_gas_is_active(j) || !bh->CaptureEnabled)
        continue;

      const double lambda = bh_ffr_capture_lambda(bh, r);
      if(!(lambda > 0))
        continue;

      if(CapturePass == BH_FFR_CAPTURE_PASS_LAMBDA)
        {
          LambdaSink[j] += lambda;
          continue;
        }

      if(!(SinkMass[j] > 0) || !(LambdaSink[j] > 0))
        continue;

      const double share = SinkMass[j] * lambda / LambdaSink[j];
      if(!isfinite(share) || share < 0 || share > SinkMass[j] * (1.0 + 1.0e-10))
        terminate("BH_FFR: invalid overlap share=%g sink=%g lambda=%g Lambda=%g", share, SinkMass[j], lambda, LambdaSink[j]);

      out.CapturedMass += share;

      const double dt_gas = bh_ffr_active_gas_timestep_code_time(j);
      if(!(dt_gas > 0) || !isfinite(dt_gas))
        terminate("BH_FFR: invalid active gas timestep=%g for supply rate ID=%llu", dt_gas,
                  (unsigned long long)P[j].ID);
      out.CapturedRate += share / dt_gas;

      double vphys2 = 0.0;
      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      for(int k = 0; k < 3; k++)
        {
          /* Transfer exactly the specific momentum represented by AREPO's
           * conserved gas variable. P[j].Vel is reconstructed from this state
           * but can differ at roundoff level after kick/update ordering; using
           * it here while removing SphP.Momentum below creates a false ledger
           * mismatch for very small sink masses. */
          const double vcons = SphP[j].Momentum[k] / P[j].Mass;
          const double dp = share * vcons;
          out.CapturedMomentum[k] += dp;
          out.CapturedAbsMomentum[k] += fabs(dp);
          const double vphys = vcons / a;
          vphys2 += vphys * vphys;
        }
      out.CapturedKineticEnergy += 0.5 * share * vphys2;

      double xtmp, ytmp, ztmp;
      double dr[3] = {NEAREST_X(P[j].Pos[0] - bh->Pos[0]), NEAREST_Y(P[j].Pos[1] - bh->Pos[1]),
                      NEAREST_Z(P[j].Pos[2] - bh->Pos[2])};
      double dv[3] = {P[j].Vel[0] - bh->Vel[0], P[j].Vel[1] - bh->Vel[1], P[j].Vel[2] - bh->Vel[2]};

      double ell[3] = {dr[1] * dv[2] - dr[2] * dv[1], dr[2] * dv[0] - dr[0] * dv[2],
                       dr[0] * dv[1] - dr[1] * dv[0]};
      const double ellnorm = sqrt(ell[0] * ell[0] + ell[1] * ell[1] + ell[2] * ell[2]);

      if(ellnorm > 0 && isfinite(ellnorm))
        for(int k = 0; k < 3; k++)
          out.CoherenceIncrement[k] += share * ell[k] / ellnorm;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static double bh_ffr_active_gas_timestep_code_time(int j)
{
  const int bin = bh_ffr_decode_hydro_timebin(P[j].TimeBinHydro);
  if(bin <= 0 || bin >= TIMEBINS || !TimeBinSynchronized[bin])
    terminate("BH_FFR: requested sink timestep for inactive/invalid gas ID=%llu bin=%d", (unsigned long long)P[j].ID, bin);

  const integertime ti_step = ((integertime)1) << bin;
  if(All.Ti_Current < ti_step)
    return 0.0;

  return bh_ffr_integer_interval_to_physical_code_time(All.Ti_Current - ti_step, All.Ti_Current);
}

static void bh_ffr_prepare_benchmark_rates(void)
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;

  for(int n = 0; n < CaptureNTargets; n++)
    {
      const int p = bh_ffr_capture_particle_from_target(n, "bh_ffr_prepare_benchmark_rates");
      const int b = P[p].BHDataIndex;
      struct bh_ffr_capture_result *res = &CaptureResults[n];

      struct bh_benchmark_environment env;
      memset(&env, 0, sizeof(env));
      env.GasMass = res->EnvGasMass;
      env.FFRRawRate = res->EnvFFRRawRate;

      if(res->EnvGasMass > 0)
        {
          if(!(res->EnvVolume > 0) || !isfinite(res->EnvVolume))
            terminate("BH_BENCHMARK: non-positive common aperture volume=%g for ID=%llu",
                      res->EnvVolume, (unsigned long long)P[p].ID);

          const double physical_volume = res->EnvVolume * a * a * a;
          env.Density = res->EnvGasMass / physical_volume;
          env.SoundSpeed = res->EnvSoundVolumeWeighted / res->EnvVolume;

          double vrel2 = 0.0;
          double j2 = 0.0;
          for(int k = 0; k < 3; k++)
            {
              const double vbulk = res->EnvVelocityVolumeWeighted[k] / res->EnvVolume;
              const double vbh = P[p].Vel[k] / a;
              const double dv = vbulk - vbh;
              vrel2 += dv * dv;
              j2 += res->EnvAngularMomentum[k] * res->EnvAngularMomentum[k];
            }
          env.RelativeSpeed = sqrt(vrel2);
          env.Vphi = sqrt(j2) / (res->EnvGasMass * All.BHAccretionRadius);

          const double rho_cgs =
              env.Density * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
          env.HydrogenNumberDensity = HYDROGEN_MASSFRAC * rho_cgs / PROTONMASS;
        }

      struct bh_benchmark_rate_result rate;
      bh_benchmark_compute_accretion(&env, BHP[b].BHMass, &rate);

      res->ModelRawRate = rate.RawRate;
      res->ModelOperationalRate = rate.OperationalRate;
      res->ModelEddingtonRate = rate.EddingtonRate;
      res->ModelBoostFactor = rate.BoostFactor;
      res->ModelAMLimiter = rate.AngularMomentumLimiter;
      res->LambdaScale = 0.0;
      res->UniformLambda = 0.0;

      if(All.BHBenchmarkAccretionModel == BH_BENCHMARK_ACC_FFR)
        {
          if(rate.RawRate > 0)
            res->LambdaScale = rate.OperationalRate / rate.RawRate;
        }
      else if(env.GasMass > 0)
        res->UniformLambda = rate.OperationalRate / env.GasMass;

      if(!isfinite(res->LambdaScale) || res->LambdaScale < 0 || res->LambdaScale > 1.0 + 1.0e-12 ||
         !isfinite(res->UniformLambda) || res->UniformLambda < 0)
        terminate("BH_BENCHMARK: invalid sink normalization scale=%g uniform=%g for ID=%llu",
                  res->LambdaScale, res->UniformLambda, (unsigned long long)P[p].ID);

      printf("BH_BENCHMARK: accretion ID=%llu task=%d model=%s target=%s "
             "Mgas=%g rho=%g cs=%g vrel=%g Vphi=%g nH=%g raw=%g edd=%g operational=%g "
             "boost=%g amlim=%g\n",
             (unsigned long long)P[p].ID, ThisTask,
             bh_benchmark_accretion_model_name(All.BHBenchmarkAccretionModel),
             All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_DIRECT ? "direct" : "reservoir",
             env.GasMass, env.Density, env.SoundSpeed, env.RelativeSpeed, env.Vphi,
             env.HydrogenNumberDensity, rate.RawRate, rate.EddingtonRate, rate.OperationalRate,
             rate.BoostFactor, rate.AngularMomentumLimiter);
      fflush(stdout);
    }
}

static void bh_ffr_compute_exact_cell_sinks(void)
{
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      const int j = TimeBinsHydro.ActiveParticleList[idx];
      if(j < 0)
        continue;
      if(j >= NumGas || P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      const double lambda = LambdaSink[j];
      if(!(lambda > 0))
        continue;

      const double dt = bh_ffr_active_gas_timestep_code_time(j);
      if(!(dt > 0))
        continue;

      const double x = lambda * dt;
      double frac = (x > 700.0) ? 1.0 : -expm1(-x);
      if(frac > All.BHMaxSinkFraction)
        frac = All.BHMaxSinkFraction;

      if(!isfinite(frac) || frac < 0 || frac >= 1)
        terminate("BH_FFR: invalid exact sink fraction=%g lambda=%g dt=%g", frac, lambda, dt);

      SinkMass[j] = P[j].Mass * frac;
    }
}

static void bh_ffr_apply_gas_sink(double *local_removed_mass, double local_removed_momentum[3],
                                  double local_removed_abs_momentum[3])
{
  *local_removed_mass = 0;
  for(int k = 0; k < 3; k++)
    {
      local_removed_momentum[k] = 0;
      local_removed_abs_momentum[k] = 0;
    }

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      const int j = TimeBinsHydro.ActiveParticleList[idx];
      if(j < 0 || j >= NumGas || !(SinkMass[j] > 0))
        continue;

      const double oldmass = P[j].Mass;
      const double dm = SinkMass[j];
      if(!(oldmass > 0) || dm < 0 || dm >= oldmass)
        terminate("BH_FFR: unsafe gas sink ID=%llu oldmass=%g dm=%g", (unsigned long long)P[j].ID, oldmass, dm);

      /* Form the removed fraction directly from dm/oldmass.  For the
       * benchmark smoke tests dm/oldmass can be extremely small; computing
       * it indirectly as 1-(oldmass-dm)/oldmass loses significant digits
       * and can create a false momentum-ledger mismatch at ~1e-9 relative. */
      const double removed_frac = dm / oldmass;
      const double keep = 1.0 - removed_frac;
      if(!(removed_frac > 0) || removed_frac >= 1 || !(keep > 0) || keep > 1)
        terminate("BH_FFR: invalid sink fractions removed=%g retained=%g for ID=%llu",
                  removed_frac, keep, (unsigned long long)P[j].ID);

      *local_removed_mass += dm;
      for(int k = 0; k < 3; k++)
        {
          const double removed_p = removed_frac * SphP[j].Momentum[k];
          local_removed_momentum[k] += removed_p;
          local_removed_abs_momentum[k] += fabs(removed_p);
          SphP[j].Momentum[k] -= removed_p;
        }

      P[j].Mass = oldmass - dm;
      SphP[j].Energy *= keep;

#ifdef PASSIVE_SCALARS
      for(int k = 0; k < PASSIVE_SCALARS; k++)
        SphP[j].PConservedScalars[k] *= keep;
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
      SphP[j].HighResMass *= keep;
#endif

      SphP[j].Density = P[j].Mass / SphP[j].Volume;
      SphP[j].OldMass = P[j].Mass;
      set_pressure_of_cell(j);

      if(!isfinite(P[j].Mass) || !(P[j].Mass > 0) || !isfinite(SphP[j].Density) || !(SphP[j].Density > 0))
        terminate("BH_FFR: invalid post-sink gas state for ID=%llu", (unsigned long long)P[j].ID);
    }
}

static void bh_ffr_check_capture_ledger(double local_removed_mass, const double local_removed_momentum[3],
                                        const double local_removed_abs_momentum[3],
                                        double local_captured_mass, const double local_captured_momentum[3],
                                        const double local_captured_abs_momentum[3])
{
  double in[14] = {local_removed_mass, local_captured_mass,
                   local_removed_momentum[0], local_removed_momentum[1], local_removed_momentum[2],
                   local_captured_momentum[0], local_captured_momentum[1], local_captured_momentum[2],
                   local_removed_abs_momentum[0], local_removed_abs_momentum[1], local_removed_abs_momentum[2],
                   local_captured_abs_momentum[0], local_captured_abs_momentum[1], local_captured_abs_momentum[2]};
  double out[14];
  MPI_Allreduce(in, out, 14, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  const double mass_scale = dmax(1.0e-30, dmax(fabs(out[0]), fabs(out[1])));
  if(fabs(out[0] - out[1]) > 1.0e-9 * mass_scale)
    terminate("BH_FFR: capture mass ledger failed removed=%g captured=%g", out[0], out[1]);

  for(int k = 0; k < 3; k++)
    {
      const double removed = out[2 + k];
      const double captured = out[5 + k];
      const double abs_removed = out[8 + k];
      const double abs_captured = out[11 + k];

      /* A signed Cartesian component can be tiny because many positive and
       * negative cell momenta cancel.  Scale roundoff against the total
       * absolute momentum transported, not the cancellation residual. */
      const double pscale =
          dmax(1.0e-30, dmax(dmax(fabs(removed), fabs(captured)),
                             dmax(abs_removed, abs_captured)));
      const double pdiff = fabs(removed - captured);

      if(pdiff > 1.0e-9 * pscale)
        terminate("BH_FFR: capture momentum ledger failed component=%d removed=%g captured=%g "
                  "absRemoved=%g absCaptured=%g diff=%g rel=%g",
                  k, removed, captured, abs_removed, abs_captured, pdiff, pdiff / pscale);
    }
}

static void bh_ffr_apply_capture_to_bhs(double *local_captured_mass, double local_captured_momentum[3],
                                        double local_captured_abs_momentum[3])
{
  *local_captured_mass = 0;
  for(int k = 0; k < 3; k++)
    {
      local_captured_momentum[k] = 0;
      local_captured_abs_momentum[k] = 0;
    }

  for(int n = 0; n < CaptureNTargets; n++)
    {
      const int p = bh_ffr_capture_particle_from_target(n, "bh_ffr_apply_capture_to_bhs");
      const int b = P[p].BHDataIndex;
      struct bh_ffr_capture_result *res = &CaptureResults[n];

      bh_ffr_set_min_neighbour_timebin(p, res->MinHydroTimeBin);

      const double dm = res->CapturedMass;
      if(!isfinite(dm) || dm < 0)
        terminate("BH_FFR: invalid captured mass=%g for particle ID=%llu", dm, (unsigned long long)P[p].ID);

      const double expected_dyn_mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
      const double dyn_scale = dmax(1.0e-30, dmax(fabs(expected_dyn_mass), fabs(P[p].Mass)));
      if(fabs(P[p].Mass - expected_dyn_mass) > 1.0e-10 * dyn_scale)
        terminate("BH_FFR: pre-capture dynamical-mass mismatch ID=%llu P.Mass=%g components=%g", (unsigned long long)P[p].ID,
                  P[p].Mass, expected_dyn_mass);

      const double old_dyn_mass = P[p].Mass;
      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
      double old_momentum[3];
      double old_vphys2 = 0.0;
      for(int k = 0; k < 3; k++)
        {
          old_momentum[k] = old_dyn_mass * P[p].Vel[k];
          const double vphys = P[p].Vel[k] / a;
          old_vphys2 += vphys * vphys;
        }
      const double old_bh_kinetic = 0.5 * old_dyn_mass * old_vphys2;

      if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_RESERVOIR)
        {
          BHP[b].ReservoirMass += dm;
          for(int k = 0; k < 3; k++)
            BHP[b].Coherence[k] += res->CoherenceIncrement[k];
        }
      else
        {
          if(BHP[b].ReservoirMass != 0 || BHP[b].WindMassBuffer != 0 ||
             BHP[b].WindMomentumBuffer != 0 || BHP[b].WindEnergyBuffer != 0 ||
             BHP[b].JetEnergyBuffer != 0)
            terminate("BH_BENCHMARK: direct accretion requires empty FFR reservoir/feedback buffers for ID=%llu",
                      (unsigned long long)P[p].ID);
          BHP[b].BHMass += dm;
        }

      P[p].Mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
      if(!(P[p].Mass > 0) || !isfinite(P[p].Mass))
        terminate("BH_FFR: invalid post-capture dynamical mass for ID=%llu", (unsigned long long)P[p].ID);

      double new_vphys2 = 0.0;
      for(int k = 0; k < 3; k++)
        {
          P[p].Vel[k] = (old_momentum[k] + res->CapturedMomentum[k]) / P[p].Mass;
          const double vphys = P[p].Vel[k] / a;
          new_vphys2 += vphys * vphys;
        }

      const double new_bh_kinetic = 0.5 * P[p].Mass * new_vphys2;
      double capture_dissipation = old_bh_kinetic + res->CapturedKineticEnergy - new_bh_kinetic;
      const double kinetic_scale =
          fmax(fabs(old_bh_kinetic) + fabs(res->CapturedKineticEnergy) + fabs(new_bh_kinetic), 1.0e-30);
      if(capture_dissipation < 0 && capture_dissipation > -5.0e-8 * kinetic_scale)
        capture_dissipation = 0.0;
      if(!isfinite(capture_dissipation) || capture_dissipation < 0)
        terminate("BH_FFR: negative/non-finite unresolved capture dissipation ID=%llu Ediss=%g scale=%g",
                  (unsigned long long)P[p].ID, capture_dissipation, kinetic_scale);

      if(dm > 0)
        {
          const double eerg = capture_dissipation * All.UnitEnergy_in_cgs / All.HubbleParam;
          printf("BH_FFR: capture dissipation ID=%llu task=%d dM=%g Ediss=%g erg\n",
                 (unsigned long long)P[p].ID, ThisTask, dm, eerg);
          fflush(stdout);
        }

      if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_RESERVOIR)
        {
          double cnorm2 = 0;
          for(int k = 0; k < 3; k++)
            cnorm2 += BHP[b].Coherence[k] * BHP[b].Coherence[k];
          const double cnorm = sqrt(cnorm2);
          const double coherence = (BHP[b].ReservoirMass > 0) ? cnorm / BHP[b].ReservoirMass : 0;

          if(!isfinite(coherence) || coherence > 1.0 + 1.0e-10)
            terminate("BH_FFR: coherence invariant failed ID=%llu C=%g Md=%g |Cd|=%g",
                      (unsigned long long)P[p].ID, coherence, BHP[b].ReservoirMass, cnorm);

          if(coherence >= All.BHMinCoherence && cnorm > 0)
            for(int k = 0; k < 3; k++)
              BHP[b].DiscDir[k] = BHP[b].Coherence[k] / cnorm;
        }

      if(!isfinite(res->CapturedRate) || res->CapturedRate < 0)
        terminate("BH_FFR: invalid captured supply rate=%g for particle ID=%llu", res->CapturedRate,
                  (unsigned long long)P[p].ID);

      /* Preserve the algebraic estimator separately from the realized
       * conservative gas sink.  TNG selects its feedback state from the
       * uncapped estimator/Eddington ratio, while its energy budget below is
       * tied to the mass that was actually captured. */
      BHP[b].BenchmarkMdotRaw = res->ModelRawRate;
      BHP[b].BenchmarkMdotOperational = res->ModelOperationalRate;

      /* Each gas sink is integrated over that gas cell's own synchronized
       * hydro step.  Summing dm_i/dt_i avoids spuriously dividing a gas-step
       * capture event by a much shorter feedback-limited BH timestep. */
      BHP[b].MdotSupply = res->CapturedRate;

      if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_DIRECT)
        {
          BHP[b].MdotProcessed = res->CapturedRate;
          BHP[b].MdotHorizon = res->CapturedRate;
          BHP[b].MdotWind = 0.0;
          BHP[b].MdotEddington = res->ModelEddingtonRate;
          BHP[b].ProcessedEddingtonRatio =
              (res->ModelEddingtonRate > 0) ? res->CapturedRate / res->ModelEddingtonRate : 0.0;
          BHP[b].ColdBlendWeight = 0.0;
          BHP[b].BolometricLuminosity = 0.0;
          BHP[b].WindPower = 0.0;
          BHP[b].JetPower = 0.0;
          BHP[b].AccretionState = BH_FFR_STATE_UNINITIALIZED;
        }

      *local_captured_mass += dm;
      for(int k = 0; k < 3; k++)
        {
          local_captured_momentum[k] += res->CapturedMomentum[k];
          local_captured_abs_momentum[k] += res->CapturedAbsMomentum[k];
        }
    }
}

void bh_ffr_capture_self_test(void)
{
  const double lambda1 = 0.3;
  const double lambda2 = 0.7;
  const double dt = 0.4;
  const double mass = 2.0;
  const double lambda = lambda1 + lambda2;
  const double sink = mass * (1.0 - exp(-lambda * dt));
  const double share1 = sink * lambda1 / lambda;
  const double share2 = sink * lambda2 / lambda;

  if(fabs((share1 + share2) - sink) > 1.0e-14 * mass)
    terminate("BH_FFR: capture self-test failed overlap partition");

  const double reversed1 = sink * lambda2 / lambda;
  const double reversed2 = sink * lambda1 / lambda;
  if(fabs((reversed1 + reversed2) - sink) > 1.0e-14 * mass)
    terminate("BH_FFR: capture self-test failed order independence");

  const double capped = dmin(sink, 0.25 * mass);
  if(!(capped > 0) || capped > 0.25 * mass)
    terminate("BH_FFR: capture self-test failed sink cap");

  /* Source-fidelity/full-unresolved-mass switch regression. For alpha=0 the
   * capture coefficient scales as lambda proportional to sqrt(M_cen), so
   * changing M_cen from M_BH to 4 M_BH must double lambda. */
  data_in bh_test;
  memset(&bh_test, 0, sizeof(bh_test));
  bh_test.CaptureEnabled = 1;
  bh_test.SchwarzschildRadius = 1.0;
  bh_test.CentralMass = 1.0;
  const double lambda_bh_only = bh_ffr_capture_lambda_core(&bh_test, 0.5, 1.0, 0.0);
  bh_test.CentralMass = 4.0;
  const double lambda_full = bh_ffr_capture_lambda_core(&bh_test, 0.5, 1.0, 0.0);

  if(!(lambda_bh_only > 0) || !isfinite(lambda_full) ||
     fabs(lambda_full / lambda_bh_only - 2.0) > 2.0e-13)
    terminate("BH_FFR: capture self-test failed central-mass switch scaling lambdaBH=%g lambdaFull=%g",
              lambda_bh_only, lambda_full);
}

void bh_ffr_capture_resolved_gas(void)
{
  /* Collective routine: ranks with zero local BH targets must still enter the
   * generic communication passes so they can receive/export neighbour-tree
   * work for BHs owned by other ranks. */
  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

#ifdef MHD
  terminate("BH_FFR: resolved gas capture is not yet enabled with MHD; a magnetic-flux sink policy is required first");
#endif

  if(!CaptureSelfTestDone)
    {
      bh_ffr_capture_self_test();
      bh_benchmark_accretion_self_test();
      CaptureSelfTestDone = 1;
    }

  CaptureNTargets = NumActiveBHFFR;
  CaptureResults = (struct bh_ffr_capture_result *)mymalloc(
      "BHFFRCaptureResults", (CaptureNTargets > 0 ? CaptureNTargets : 1) * sizeof(*CaptureResults));
  if(CaptureNTargets > 0)
    memset(CaptureResults, 0, CaptureNTargets * sizeof(*CaptureResults));
  for(int n = 0; n < CaptureNTargets; n++)
    CaptureResults[n].MinHydroTimeBin = TIMEBINS;

  LambdaSink = NULL;
  SinkMass = NULL;
  if(NumGas > 0)
    {
      LambdaSink = (double *)mymalloc("BHFFRLambdaSink", NumGas * sizeof(double));
      SinkMass = (double *)mymalloc("BHFFRSinkMass", NumGas * sizeof(double));
      memset(LambdaSink, 0, NumGas * sizeof(double));
      memset(SinkMass, 0, NumGas * sizeof(double));
    }

  if(All.TotNumGas > 0)
    {
      CapturePass = BH_FFR_CAPTURE_PASS_ENVIRONMENT;
      generic_set_MaxNexport();
      generic_comm_pattern(CaptureNTargets, kernel_local, kernel_imported);

      bh_ffr_prepare_benchmark_rates();

      CapturePass = BH_FFR_CAPTURE_PASS_LAMBDA;
      generic_set_MaxNexport();
      generic_comm_pattern(CaptureNTargets, kernel_local, kernel_imported);

      bh_ffr_compute_exact_cell_sinks();

      CapturePass = BH_FFR_CAPTURE_PASS_SHARES;
      generic_set_MaxNexport();
      generic_comm_pattern(CaptureNTargets, kernel_local, kernel_imported);
    }

  for(int n = 0; n < CaptureNTargets; n++)
    if(CaptureResults[n].MinHydroTimeBin == TIMEBINS)
      CaptureResults[n].MinHydroTimeBin = -1;

  double local_captured_mass, local_captured_momentum[3], local_captured_abs_momentum[3];
  bh_ffr_apply_capture_to_bhs(&local_captured_mass, local_captured_momentum, local_captured_abs_momentum);

  double local_removed_mass, local_removed_momentum[3], local_removed_abs_momentum[3];
  bh_ffr_apply_gas_sink(&local_removed_mass, local_removed_momentum, local_removed_abs_momentum);

  bh_ffr_check_capture_ledger(local_removed_mass, local_removed_momentum, local_removed_abs_momentum,
                              local_captured_mass, local_captured_momentum, local_captured_abs_momentum);

  if(SinkMass != NULL)
    myfree(SinkMass);
  if(LambdaSink != NULL)
    myfree(LambdaSink);
  myfree(CaptureResults);

  SinkMass = NULL;
  LambdaSink = NULL;
  CaptureResults = NULL;
  CaptureNTargets = 0;
}
