#include <float.h>
#include <math.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/* Exact finite-step reservoir processing plus hysteretic
 * ADIOS/truncated/cold state classification.  Iteration 6 now commits the
 * candidate only after bh_ffr_apply_inner_flow() has partitioned it exactly
 * into horizon and wind channels. */

static double bh_ffr_eddington_rate_cgs(double mass_g, double efficiency)
{
  if(!(mass_g > 0))
    return 0.0;
  if(!(efficiency > 0) || !(efficiency < 1))
    terminate("BH_FFR: invalid Eddington reference efficiency=%g", efficiency);

  const double rate = 4.0 * M_PI * GRAVITY * mass_g * PROTONMASS / (efficiency * THOMPSON * CLIGHT);
  if(!isfinite(rate) || !(rate > 0))
    terminate("BH_FFR: invalid cgs Eddington rate=%g for mass=%g g", rate, mass_g);

  return rate;
}

double bh_ffr_reservoir_processable_mass(double reservoir_mass, double dt_myr, double tau_myr)
{
  if(!isfinite(reservoir_mass) || reservoir_mass < 0)
    terminate("BH_FFR: invalid reservoir mass=%g in processable-mass calculation", reservoir_mass);
  if(!isfinite(dt_myr) || dt_myr < 0)
    terminate("BH_FFR: invalid reservoir timestep=%g Myr", dt_myr);
  if(!isfinite(tau_myr) || !(tau_myr > 0))
    terminate("BH_FFR: invalid reservoir processing time=%g Myr", tau_myr);

  if(reservoir_mass == 0 || dt_myr == 0)
    return 0.0;

  /* -expm1(-x) is accurate for x << 1 and tends safely to unity for x >> 1. */
  double processed = reservoir_mass * (-expm1(-dt_myr / tau_myr));
  if(processed < 0 && processed > -1.0e-14 * dmax(reservoir_mass, DBL_MIN))
    processed = 0;
  if(processed > reservoir_mass && processed < reservoir_mass * (1.0 + 1.0e-14))
    processed = reservoir_mass;

  if(!isfinite(processed) || processed < 0 || processed > reservoir_mass)
    terminate("BH_FFR: invalid processable reservoir mass=%g from Md=%g dt=%g tau=%g", processed, reservoir_mass, dt_myr,
              tau_myr);

  return processed;
}

double bh_ffr_eddington_rate_code(double bh_mass)
{
  if(!isfinite(bh_mass) || bh_mass < 0)
    terminate("BH_FFR: invalid BH mass=%g in Eddington-rate calculation", bh_mass);
  if(bh_mass == 0)
    return 0.0;

  if(!(All.UnitMass_in_g > 0) || !(All.UnitTime_in_s > 0) || !(All.HubbleParam > 0))
    terminate("BH_FFR: invalid units for Eddington conversion Munit=%g Tunit=%g h=%g", All.UnitMass_in_g, All.UnitTime_in_s,
              All.HubbleParam);

  /* Stored cosmological masses are in code mass units with the usual h
   * scaling. Convert to physical grams for L_Edd, then convert g/s back to
   * code mass per physical code-time. The latter has no residual h factor,
   * matching MdotSupply = DeltaM_code / Delta t_physical,code. */
  const double mass_g = bh_mass * All.UnitMass_in_g / All.HubbleParam;
  const double rate_cgs = bh_ffr_eddington_rate_cgs(mass_g, All.BHRadiativeEfficiency);
  const double rate_code = rate_cgs * All.UnitTime_in_s / All.UnitMass_in_g;

  if(!isfinite(rate_code) || !(rate_code > 0))
    terminate("BH_FFR: invalid internal Eddington rate=%g for BH mass=%g", rate_code, bh_mass);

  return rate_code;
}

double bh_ffr_truncation_to_hot_radius_ratio(double processed_edd_ratio)
{
  if(!isfinite(processed_edd_ratio) || processed_edd_ratio < 0)
    terminate("BH_FFR: invalid processed Eddington ratio=%g", processed_edd_ratio);
  if(!(All.BHRHotMaxInRs > 3.0))
    terminate("BH_FFR: BHRHotMaxInRs=%g must exceed 3", All.BHRHotMaxInRs);

  if(processed_edd_ratio == 0)
    return DBL_MAX;

  const double prefactor = 3.0 / All.BHRHotMaxInRs;
  const double q = BH_FFR_COLD_NOMINAL_EDD_RATIO / processed_edd_ratio;

  /* Avoid overflowing DBL_MAX/prefactor when prefactor < 1. */
  if(q > sqrt(DBL_MAX) / sqrt(prefactor))
    return DBL_MAX;

  const double ratio = prefactor * q * q;
  if(!isfinite(ratio) || ratio < 0)
    terminate("BH_FFR: invalid Rtr/Rhot=%g at processed Eddington ratio=%g", ratio, processed_edd_ratio);

  return ratio;
}

static int bh_ffr_nominal_accretion_state(double processed_edd_ratio, double rtr_over_rhot)
{
  if(processed_edd_ratio >= BH_FFR_COLD_NOMINAL_EDD_RATIO)
    return BH_FFR_STATE_COLD;
  if(rtr_over_rhot >= 1.0)
    return BH_FFR_STATE_ADIOS;
  return BH_FFR_STATE_TRUNCATED;
}

int bh_ffr_classify_accretion_state(int previous_state, double processed_edd_ratio, double rtr_over_rhot)
{
  if(!isfinite(processed_edd_ratio) || processed_edd_ratio < 0 || !isfinite(rtr_over_rhot) || rtr_over_rhot < 0)
    terminate("BH_FFR: invalid state-classification inputs dotm=%g Rtr/Rhot=%g", processed_edd_ratio, rtr_over_rhot);

  if(previous_state == BH_FFR_STATE_UNINITIALIZED)
    return bh_ffr_nominal_accretion_state(processed_edd_ratio, rtr_over_rhot);

  if(previous_state < BH_FFR_STATE_ADIOS || previous_state > BH_FFR_STATE_COLD)
    terminate("BH_FFR: invalid previous accretion state=%d", previous_state);

  if(previous_state == BH_FFR_STATE_COLD)
    {
      if(processed_edd_ratio >= BH_FFR_COLD_EXIT_EDD_RATIO)
        return BH_FFR_STATE_COLD;
      return (rtr_over_rhot > BH_FFR_ADIOS_ENTER_RTR_FACTOR) ? BH_FFR_STATE_ADIOS : BH_FFR_STATE_TRUNCATED;
    }

  if(previous_state == BH_FFR_STATE_ADIOS)
    {
      if(rtr_over_rhot >= BH_FFR_ADIOS_EXIT_RTR_FACTOR)
        return BH_FFR_STATE_ADIOS;
      if(processed_edd_ratio >= BH_FFR_COLD_ENTER_EDD_RATIO)
        return BH_FFR_STATE_COLD;
      return BH_FFR_STATE_TRUNCATED;
    }

  if(processed_edd_ratio >= BH_FFR_COLD_ENTER_EDD_RATIO)
    return BH_FFR_STATE_COLD;
  if(rtr_over_rhot > BH_FFR_ADIOS_ENTER_RTR_FACTOR)
    return BH_FFR_STATE_ADIOS;
  return BH_FFR_STATE_TRUNCATED;
}

void bh_ffr_update_reservoir_state(int p, double dt_myr, double dt_code)
{
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid particle index=%d in reservoir-state update", p);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for particle ID=%llu in reservoir-state update", (unsigned long long)P[p].ID);

  if(!isfinite(dt_myr) || dt_myr < 0 || !isfinite(dt_code) || dt_code < 0)
    terminate("BH_FFR: invalid reservoir-state timestep dt_myr=%g dt_code=%g", dt_myr, dt_code);
  if((dt_myr == 0) != (dt_code == 0))
    terminate("BH_FFR: inconsistent zero timestep dt_myr=%g dt_code=%g", dt_myr, dt_code);

  const double processable = bh_ffr_reservoir_processable_mass(BHP[b].ReservoirMass, dt_myr, All.BHDiskTimeMyr);
  const double mdot_processed = (dt_code > 0) ? processable / dt_code : 0.0;
  const double mdot_edd = bh_ffr_eddington_rate_code(BHP[b].BHMass);
  const double processed_edd_ratio = (mdot_edd > 0) ? mdot_processed / mdot_edd : 0.0;
  const double rtr_over_rhot = bh_ffr_truncation_to_hot_radius_ratio(processed_edd_ratio);
  const int state = bh_ffr_classify_accretion_state(BHP[b].AccretionState, processed_edd_ratio, rtr_over_rhot);

  if(!isfinite(mdot_processed) || mdot_processed < 0 || !isfinite(mdot_edd) || mdot_edd < 0 ||
     !isfinite(processed_edd_ratio) || processed_edd_ratio < 0)
    terminate("BH_FFR: invalid reservoir/inner-flow diagnostics for particle ID=%llu", (unsigned long long)P[p].ID);

  BHP[b].MdotProcessed = mdot_processed;
  BHP[b].MdotEddington = mdot_edd;
  BHP[b].ProcessedEddingtonRatio = processed_edd_ratio;
  BHP[b].AccretionState = state;

  /* The transaction is atomic at the sub-grid level: no reservoir mass is
   * removed unless the same call assigns it to the horizon or wind channel. */
  bh_ffr_apply_inner_flow(p, processable, dt_code);
}

void bh_ffr_reservoir_self_test(void)
{
  const double md = 2.0;
  const double exact = md * (1.0 - exp(-1.0));
  const double got = bh_ffr_reservoir_processable_mass(md, 5.0, 5.0);
  if(fabs(got - exact) > 1.0e-14 * md)
    terminate("BH_FFR: reservoir self-test failed exact exponential got=%g expected=%g", got, exact);

  const double tiny = bh_ffr_reservoir_processable_mass(md, 1.0e-9, 5.0);
  const double tiny_linear = md * 1.0e-9 / 5.0;
  if(fabs(tiny - tiny_linear) > 1.0e-9 * tiny_linear)
    terminate("BH_FFR: reservoir self-test failed small-step limit got=%g expected=%g", tiny, tiny_linear);

  if(bh_ffr_reservoir_processable_mass(0.0, 5.0, 5.0) != 0.0 ||
     bh_ffr_reservoir_processable_mass(md, 0.0, 5.0) != 0.0)
    terminate("BH_FFR: reservoir self-test failed zero-mass/zero-step limit");

  const double saturated = bh_ffr_reservoir_processable_mass(md, 1.0e6, 5.0);
  if(fabs(saturated - md) > 1.0e-14 * md)
    terminate("BH_FFR: reservoir self-test failed large-step saturation got=%g expected=%g", saturated, md);

  const double edd_cgs = bh_ffr_eddington_rate_cgs(1.0e5 * SOLAR_MASS, 0.1);
  const double edd_msun_yr = edd_cgs * SEC_PER_YEAR / SOLAR_MASS;
  if(fabs(edd_msun_yr - 2.21963658503e-3) > 1.0e-8 * 2.21963658503e-3)
    terminate("BH_FFR: reservoir self-test failed Eddington normalization got=%g Msun/yr", edd_msun_yr);

  const double edd_cgs_twice = bh_ffr_eddington_rate_cgs(2.0e5 * SOLAR_MASS, 0.1);
  if(fabs(edd_cgs_twice - 2.0 * edd_cgs) > 1.0e-14 * edd_cgs_twice)
    terminate("BH_FFR: reservoir self-test failed linear Eddington mass scaling");

  const double boundary_ratio =
      (3.0 / 1200.0) * pow(BH_FFR_COLD_NOMINAL_EDD_RATIO / 1.0e-3, 2.0);
  if(fabs(boundary_ratio - 1.0) > 1.0e-14)
    terminate("BH_FFR: reservoir self-test failed ADIOS boundary identity Rtr/Rhot=%g", boundary_ratio);

  if(bh_ffr_classify_accretion_state(BH_FFR_STATE_UNINITIALIZED, 0.0, DBL_MAX) != BH_FFR_STATE_ADIOS ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_UNINITIALIZED, 0.01, 0.01) != BH_FFR_STATE_TRUNCATED ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_UNINITIALIZED, 0.03, 0.001) != BH_FFR_STATE_COLD)
    terminate("BH_FFR: reservoir self-test failed nominal state classification");

  if(bh_ffr_classify_accretion_state(BH_FFR_STATE_TRUNCATED, 0.021, 0.0025) != BH_FFR_STATE_TRUNCATED ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_TRUNCATED, 0.023, 0.0025) != BH_FFR_STATE_COLD ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_COLD, 0.019, 0.003) != BH_FFR_STATE_COLD ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_COLD, 0.017, 0.004) != BH_FFR_STATE_TRUNCATED ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_COLD, 9.0e-4, 1.21) != BH_FFR_STATE_ADIOS ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_COLD, 9.5e-4, 1.10) != BH_FFR_STATE_TRUNCATED)
    terminate("BH_FFR: reservoir self-test failed cold-state hysteresis");

  if(bh_ffr_classify_accretion_state(BH_FFR_STATE_TRUNCATED, 9.0e-4, 1.21) != BH_FFR_STATE_ADIOS ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_ADIOS, 1.05e-3, 0.90) != BH_FFR_STATE_ADIOS ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_ADIOS, 1.13e-3, 0.79) != BH_FFR_STATE_TRUNCATED)
    terminate("BH_FFR: reservoir self-test failed ADIOS hysteresis");

  /* Equality stays on the retained side of each hysteresis band. */
  if(bh_ffr_classify_accretion_state(BH_FFR_STATE_TRUNCATED, BH_FFR_COLD_ENTER_EDD_RATIO, 0.01) != BH_FFR_STATE_COLD ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_COLD, BH_FFR_COLD_EXIT_EDD_RATIO, 0.01) != BH_FFR_STATE_COLD ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_TRUNCATED, 9.0e-4, BH_FFR_ADIOS_ENTER_RTR_FACTOR) !=
         BH_FFR_STATE_TRUNCATED ||
     bh_ffr_classify_accretion_state(BH_FFR_STATE_ADIOS, 1.1e-3, BH_FFR_ADIOS_EXIT_RTR_FACTOR) != BH_FFR_STATE_ADIOS)
    terminate("BH_FFR: reservoir self-test failed hysteresis equality convention");
}
