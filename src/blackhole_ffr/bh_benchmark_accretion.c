#include <float.h>
#include <math.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Accretion-only comparison layer used by the convergence benchmark branch.
 *
 * The environment is deliberately supplied by the common fixed proper
 * aperture in bh_ffr_capture.c so Bondi-family models and FFR see the same
 * resolved gas.  This file contains only algebraic model closures; it does
 * not perform gas removal or feedback.
 */

static double bh_benchmark_tng_bondi_core(double g, double mass, double rho, double cs)
{
  if(!(mass > 0) || !(rho > 0))
    return 0.0;
  if(!(g > 0) || !(cs > 0) || !isfinite(g) || !isfinite(mass) || !isfinite(rho) || !isfinite(cs))
    terminate("BH_BENCHMARK: invalid TNG-Bondi inputs G=%g M=%g rho=%g cs=%g", g, mass, rho, cs);

  const double rate = 4.0 * M_PI * g * g * mass * mass * rho / (cs * cs * cs);
  if(!isfinite(rate) || rate < 0)
    terminate("BH_BENCHMARK: invalid TNG-Bondi rate=%g", rate);
  return rate;
}

static double bh_benchmark_bhl_core(double g, double mass, double rho, double cs, double vrel)
{
  if(!(mass > 0) || !(rho > 0))
    return 0.0;
  if(!(g > 0) || cs < 0 || vrel < 0 || !isfinite(g) || !isfinite(mass) || !isfinite(rho) ||
     !isfinite(cs) || !isfinite(vrel))
    terminate("BH_BENCHMARK: invalid BHL inputs G=%g M=%g rho=%g cs=%g vrel=%g", g, mass, rho, cs, vrel);

  const double q = cs * cs + vrel * vrel;
  if(!(q > 0) || !isfinite(q))
    terminate("BH_BENCHMARK: singular BHL denominator cs=%g vrel=%g", cs, vrel);

  const double denom = q * sqrt(q);
  const double rate = 4.0 * M_PI * g * g * mass * mass * rho / denom;
  if(!isfinite(rate) || rate < 0)
    terminate("BH_BENCHMARK: invalid BHL rate=%g", rate);
  return rate;
}

static double bh_benchmark_density_boost_core(double n_h, int mode, double alpha0, double nstar, double beta)
{
  if(!isfinite(n_h) || n_h < 0 || !isfinite(alpha0) || !(alpha0 > 0) ||
     !isfinite(nstar) || !(nstar > 0) || !isfinite(beta) || beta < 0)
    terminate("BH_BENCHMARK: invalid boost inputs nH=%g mode=%d alpha0=%g nstar=%g beta=%g",
              n_h, mode, alpha0, nstar, beta);

  if(mode == BH_BENCHMARK_BOOST_CONSTANT)
    return alpha0;

  if(mode == BH_BENCHMARK_BOOST_DENSITY)
    {
      if(n_h < nstar)
        return 1.0;
      const double alpha = pow(n_h / nstar, beta);
      if(!isfinite(alpha) || alpha < 1.0)
        terminate("BH_BENCHMARK: invalid density boost alpha=%g", alpha);
      return alpha;
    }

  terminate("BH_BENCHMARK: unknown boost mode=%d", mode);
  return 1.0;
}

static double bh_benchmark_am_limiter_core(double cs, double vphi, double cvisc)
{
  if(!isfinite(cs) || cs < 0 || !isfinite(vphi) || vphi < 0 || !isfinite(cvisc) || !(cvisc > 0))
    terminate("BH_BENCHMARK: invalid AM limiter inputs cs=%g Vphi=%g Cvisc=%g", cs, vphi, cvisc);

  if(vphi == 0)
    return 1.0;

  const double x = cs / vphi;
  const double limiter = (x > 0) ? (x * x * x / cvisc) : 0.0;
  if(!isfinite(limiter) || limiter < 0)
    terminate("BH_BENCHMARK: invalid AM limiter=%g", limiter);

  return dmin(1.0, limiter);
}

const char *bh_benchmark_accretion_model_name(int model)
{
  switch(model)
    {
      case BH_BENCHMARK_ACC_TNG_BONDI:
        return "tng-bondi";
      case BH_BENCHMARK_ACC_BOOSTED_BONDI:
        return "boosted-bondi";
      case BH_BENCHMARK_ACC_AM_BONDI:
        return "am-bondi";
      case BH_BENCHMARK_ACC_FFR:
        return "ffr";
      case BH_BENCHMARK_ACC_FFR_SHELL:
        return "ffr-shell";
      case BH_BENCHMARK_ACC_CONVJ_SHELL_FFR:
        return "convj-shell-ffr";
      default:
        return "unknown";
    }
}

double bh_benchmark_eddington_rate_code(double bh_mass)
{
  if(!isfinite(bh_mass) || bh_mass < 0)
    terminate("BH_BENCHMARK: invalid BH mass=%g in Eddington rate", bh_mass);
  if(bh_mass == 0)
    return 0.0;

  if(!(All.BHBenchmarkRadiativeEfficiency > 0) || !(All.BHBenchmarkRadiativeEfficiency < 1) ||
     !(All.UnitMass_in_g > 0) || !(All.UnitTime_in_s > 0) || !(All.HubbleParam > 0))
    terminate("BH_BENCHMARK: invalid Eddington parameters eps=%g Munit=%g Tunit=%g h=%g",
              All.BHBenchmarkRadiativeEfficiency, All.UnitMass_in_g, All.UnitTime_in_s, All.HubbleParam);

  const double mass_g = bh_mass * All.UnitMass_in_g / All.HubbleParam;
  const double rate_cgs =
      4.0 * M_PI * GRAVITY * mass_g * PROTONMASS /
      (All.BHBenchmarkRadiativeEfficiency * THOMPSON * CLIGHT);
  const double rate_code = rate_cgs * All.UnitTime_in_s / All.UnitMass_in_g;

  if(!isfinite(rate_code) || !(rate_code > 0))
    terminate("BH_BENCHMARK: invalid Eddington rate=%g for M=%g", rate_code, bh_mass);

  return rate_code;
}

void bh_benchmark_compute_all_raw_rates(const struct bh_benchmark_environment *env, double bh_mass,
                                        double raw_rates[BH_BENCHMARK_ACC_COUNT],
                                        double *boost_factor, double *am_limiter,
                                        double *env_factor)
{
  if(env == NULL || raw_rates == NULL || boost_factor == NULL || am_limiter == NULL ||
     env_factor == NULL)
    terminate("BH_BENCHMARK_ALL: null diagnostic input/output");

  for(int model = 0; model < BH_BENCHMARK_ACC_COUNT; model++)
    raw_rates[model] = 0.0;
  *boost_factor = 1.0;
  *am_limiter = 1.0;
  *env_factor = 1.0;

  if(bh_mass <= 0)
    return;

  /* FFR-family diagnostics are defined by their own apertures.  Do not erase
   * a valid selected shell merely because the historical common Bondi aperture
   * happens to contain no gas. */
  raw_rates[BH_BENCHMARK_ACC_FFR] = env->FFRRawRate;
  raw_rates[BH_BENCHMARK_ACC_FFR_SHELL] = env->FFRShellRawRate;
  raw_rates[BH_BENCHMARK_ACC_CONVJ_SHELL_FFR] = env->FFRConvJShellRate;
  *env_factor = env->ConvJFactor;

  if(env->GasMass <= 0)
    return;

  raw_rates[BH_BENCHMARK_ACC_TNG_BONDI] =
      bh_benchmark_tng_bondi_core(All.G, bh_mass, env->Density, env->SoundSpeed);

  const double bhl =
      bh_benchmark_bhl_core(All.G, bh_mass, env->Density, env->SoundSpeed, env->RelativeSpeed);

  *boost_factor =
      bh_benchmark_density_boost_core(env->HydrogenNumberDensity, All.BHBenchmarkBoostMode,
                                      All.BHBenchmarkBoostAlpha, All.BHBenchmarkBoostDensityThreshold,
                                      All.BHBenchmarkBoostBeta);
  raw_rates[BH_BENCHMARK_ACC_BOOSTED_BONDI] = (*boost_factor) * bhl;

  *am_limiter =
      bh_benchmark_am_limiter_core(env->SoundSpeed, env->Vphi, All.BHBenchmarkAMViscosity);
  raw_rates[BH_BENCHMARK_ACC_AM_BONDI] = (*am_limiter) * bhl;

  for(int model = 0; model < BH_BENCHMARK_ACC_COUNT; model++)
    if(!isfinite(raw_rates[model]) || raw_rates[model] < 0)
      terminate("BH_BENCHMARK_ALL: invalid raw rate model=%d rate=%g", model, raw_rates[model]);
}

void bh_benchmark_compute_accretion(const struct bh_benchmark_environment *env, double bh_mass,
                                    struct bh_benchmark_rate_result *out)
{
  if(env == NULL || out == NULL)
    terminate("BH_BENCHMARK: null accretion input/output");
  if(!isfinite(bh_mass) || bh_mass < 0 || !isfinite(env->GasMass) || env->GasMass < 0 ||
     !isfinite(env->Density) || env->Density < 0 || !isfinite(env->SoundSpeed) || env->SoundSpeed < 0 ||
     !isfinite(env->RelativeSpeed) || env->RelativeSpeed < 0 || !isfinite(env->Vphi) || env->Vphi < 0 ||
     !isfinite(env->HydrogenNumberDensity) || env->HydrogenNumberDensity < 0 ||
     !isfinite(env->FFRRawRate) || env->FFRRawRate < 0 ||
     !isfinite(env->FFRShellRawRate) || env->FFRShellRawRate < 0 ||
     !isfinite(env->FFRShellGeometricRate) || env->FFRShellGeometricRate < 0 ||
     !isfinite(env->FFRConvJShellRate) || env->FFRConvJShellRate < 0 ||
     !isfinite(env->ConvJFactor) || env->ConvJFactor < 0 || env->ConvJFactor > 1.0 + 1.0e-12 ||
     !isfinite(env->MachFactor) || env->MachFactor < 0 || env->MachFactor > 1.0 + 1.0e-12 ||
     !isfinite(env->AngularMomentumFactor) || env->AngularMomentumFactor < 0 ||
     env->AngularMomentumFactor > 1.0 + 1.0e-12)
    terminate("BH_BENCHMARK: invalid common environment Mgas=%g rho=%g cs=%g vrel=%g Vphi=%g nH=%g "
              "ffr=%g shell=%g shellGeom=%g convj=%g Fcorr=%g Fmach=%g Fj=%g",
              env->GasMass, env->Density, env->SoundSpeed, env->RelativeSpeed, env->Vphi,
              env->HydrogenNumberDensity, env->FFRRawRate, env->FFRShellRawRate,
              env->FFRShellGeometricRate, env->FFRConvJShellRate, env->ConvJFactor,
              env->MachFactor, env->AngularMomentumFactor);

  memset(out, 0, sizeof(*out));
  out->BoostFactor = 1.0;
  out->AngularMomentumLimiter = 1.0;
  out->ConvJFactor = env->ConvJFactor;
  out->MachFactor = env->MachFactor;
  out->AngularMomentumFactor = env->AngularMomentumFactor;
  out->EddingtonRate = bh_benchmark_eddington_rate_code(bh_mass);

  if(bh_mass <= 0)
    return;

  if(env->GasMass <= 0 &&
     All.BHBenchmarkAccretionModel != BH_BENCHMARK_ACC_FFR_SHELL &&
     All.BHBenchmarkAccretionModel != BH_BENCHMARK_ACC_CONVJ_SHELL_FFR)
    return;

  switch(All.BHBenchmarkAccretionModel)
    {
      case BH_BENCHMARK_ACC_TNG_BONDI:
        out->RawRate = bh_benchmark_tng_bondi_core(All.G, bh_mass, env->Density, env->SoundSpeed);
        break;

      case BH_BENCHMARK_ACC_BOOSTED_BONDI:
        out->BoostFactor =
            bh_benchmark_density_boost_core(env->HydrogenNumberDensity, All.BHBenchmarkBoostMode,
                                            All.BHBenchmarkBoostAlpha, All.BHBenchmarkBoostDensityThreshold,
                                            All.BHBenchmarkBoostBeta);
        out->RawRate = out->BoostFactor *
                       bh_benchmark_bhl_core(All.G, bh_mass, env->Density, env->SoundSpeed, env->RelativeSpeed);
        break;

      case BH_BENCHMARK_ACC_AM_BONDI:
        out->AngularMomentumLimiter =
            bh_benchmark_am_limiter_core(env->SoundSpeed, env->Vphi, All.BHBenchmarkAMViscosity);
        out->RawRate = out->AngularMomentumLimiter *
                       bh_benchmark_bhl_core(All.G, bh_mass, env->Density, env->SoundSpeed, env->RelativeSpeed);
        break;

      case BH_BENCHMARK_ACC_FFR:
        out->RawRate = env->FFRRawRate;
        break;

      case BH_BENCHMARK_ACC_FFR_SHELL:
        out->RawRate = env->FFRShellRawRate;
        break;

      case BH_BENCHMARK_ACC_CONVJ_SHELL_FFR:
        out->RawRate = env->FFRConvJShellRate;
        break;

      default:
        terminate("BH_BENCHMARK: unsupported accretion model=%d", All.BHBenchmarkAccretionModel);
    }

  if(!isfinite(out->RawRate) || out->RawRate < 0)
    terminate("BH_BENCHMARK: invalid raw accretion rate=%g", out->RawRate);

  out->OperationalRate = out->RawRate;
  if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_DIRECT)
    {
      const double cap = All.BHBenchmarkEddingtonFactor * out->EddingtonRate;
      if(!isfinite(cap) || cap < 0)
        terminate("BH_BENCHMARK: invalid direct Eddington cap=%g", cap);
      out->OperationalRate = dmin(out->RawRate, cap);
    }

  if(!isfinite(out->OperationalRate) || out->OperationalRate < 0 ||
     out->OperationalRate > out->RawRate * (1.0 + 1.0e-12))
    terminate("BH_BENCHMARK: invalid operational rate raw=%g operational=%g",
              out->RawRate, out->OperationalRate);
}

void bh_benchmark_accretion_self_test(void)
{
  const double g = 2.0, mass = 3.0, rho = 5.0, cs = 7.0, v = 11.0;
  const double b = bh_benchmark_tng_bondi_core(g, mass, rho, cs);
  const double expected_b = 4.0 * M_PI * g * g * mass * mass * rho / (cs * cs * cs);
  if(fabs(b / expected_b - 1.0) > 2.0e-14)
    terminate("BH_BENCHMARK: TNG-Bondi self-test failed got=%g expected=%g", b, expected_b);

  const double b4 = bh_benchmark_tng_bondi_core(g, 2.0 * mass, rho, cs);
  if(fabs(b4 / b - 4.0) > 2.0e-14)
    terminate("BH_BENCHMARK: Bondi M^2 scaling self-test failed");

  const double boost_low = bh_benchmark_density_boost_core(0.05, BH_BENCHMARK_BOOST_DENSITY, 100.0, 0.1, 2.0);
  const double boost_high = bh_benchmark_density_boost_core(1.0, BH_BENCHMARK_BOOST_DENSITY, 100.0, 0.1, 2.0);
  const double boost_const = bh_benchmark_density_boost_core(0.0, BH_BENCHMARK_BOOST_CONSTANT, 100.0, 0.1, 2.0);
  if(boost_low != 1.0 || fabs(boost_high - 100.0) > 1.0e-13 || boost_const != 100.0)
    terminate("BH_BENCHMARK: boosted-Bondi self-test failed low=%g high=%g const=%g",
              boost_low, boost_high, boost_const);

  const double l0 = bh_benchmark_am_limiter_core(cs, 0.0, 2.0 * M_PI);
  const double vphi = 2.0 * cs;
  const double l1 = bh_benchmark_am_limiter_core(cs, vphi, 2.0 * M_PI);
  const double l1_expected = 1.0 / (8.0 * 2.0 * M_PI);
  if(l0 != 1.0 || fabs(l1 / l1_expected - 1.0) > 2.0e-14)
    terminate("BH_BENCHMARK: AM-Bondi self-test failed l0=%g l1=%g expected=%g", l0, l1, l1_expected);

  const double bhl = bh_benchmark_bhl_core(g, mass, rho, cs, v);
  if(!(bhl > 0) || !(bhl < b))
    terminate("BH_BENCHMARK: BHL relative-velocity self-test failed BHL=%g TNG=%g", bhl, b);


}