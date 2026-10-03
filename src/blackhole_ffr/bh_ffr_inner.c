#include <float.h>
#include <math.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/* Iteration 6: quasi-steady inner-flow partition and scalar energetics.
 *
 * The newer FFR-MACER specification deliberately has no hard Eddington cap
 * on mass flow.  MdotProcessed is partitioned exactly into a horizon rest-mass
 * channel and an unresolved wind channel.  The cold branch uses the monotonic
 * Gofford mass-conservation solve rather than combining an independent fixed
 * retention fraction with an unrelated wind-loss fit.
 *
 * This file only generates the unresolved mechanical buffers. Iteration 7
 * may consume wind mass/momentum/energy in resolved bipolar packets; jet
 * energy remains buffered for the subsequent jet-feedback stage. */

#define BH_FFR_GOFFORD_MDOT_NORM_MSUN_YR 0.28
#define BH_FFR_GOFFORD_MDOT_EXPONENT 0.85
#define BH_FFR_GOFFORD_V_NORM_KMS 2.5e4
#define BH_FFR_GOFFORD_V_EXPONENT 0.4
#define BH_FFR_GOFFORD_V_MAX_KMS 1.0e5
#define BH_FFR_GOFFORD_L_NORM_CGS 1.0e45
#define BH_FFR_HOT_WIND_VK_FACTOR 0.2
#define BH_FFR_HOT_JET_EFFICIENCY 0.125

static double bh_ffr_rate_cgs_to_code(double rate_cgs)
{
  if(!isfinite(rate_cgs) || rate_cgs < 0)
    terminate("BH_FFR: invalid cgs mass rate=%g", rate_cgs);
  if(!(All.UnitMass_in_g > 0) || !(All.UnitTime_in_s > 0))
    terminate("BH_FFR: invalid mass/time units M=%g T=%g", All.UnitMass_in_g, All.UnitTime_in_s);

  return rate_cgs * All.UnitTime_in_s / All.UnitMass_in_g;
}

static double bh_ffr_power_code_to_cgs(double power_code)
{
  if(!isfinite(power_code) || power_code < 0)
    terminate("BH_FFR: invalid code power=%g", power_code);
  if(!(All.UnitEnergy_in_cgs > 0) || !(All.UnitTime_in_s > 0))
    terminate("BH_FFR: invalid energy/time units E=%g T=%g", All.UnitEnergy_in_cgs, All.UnitTime_in_s);

  return power_code * All.UnitEnergy_in_cgs / All.UnitTime_in_s;
}

static double bh_ffr_radiative_luminosity_code(double mdot_h, double mdot_edd)
{
  if(!isfinite(mdot_h) || mdot_h < 0 || !isfinite(mdot_edd) || mdot_edd < 0)
    terminate("BH_FFR: invalid horizon/Eddington rates H=%g Edd=%g in luminosity", mdot_h, mdot_edd);
  if(mdot_h == 0)
    return 0.0;
  if(!(mdot_edd > 0))
    terminate("BH_FFR: positive horizon rate=%g with non-positive Eddington rate=%g", mdot_h, mdot_edd);
  if(!(All.BHRadiativeEfficiency > 0) || !(All.BHRadiativeEfficiency < 1))
    terminate("BH_FFR: invalid radiative efficiency=%g", All.BHRadiativeEfficiency);

  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  const double dotm_h = mdot_h / mdot_edd;
  double luminosity;

  if(dotm_h <= 1.0)
    luminosity = All.BHRadiativeEfficiency * mdot_h * c_internal * c_internal * fmin(1.0, 10.0 * dotm_h);
  else
    {
      const double l_edd = All.BHRadiativeEfficiency * mdot_edd * c_internal * c_internal;
      luminosity = l_edd * (1.0 + log(dotm_h));
    }

  if(!isfinite(luminosity) || luminosity < 0)
    terminate("BH_FFR: invalid bolometric luminosity=%g at dotm_H=%g", luminosity, dotm_h);

  return luminosity;
}

/* Efficient cold-disc luminosity used only inside the Gofford mass-loss
 * closure.  The final reported luminosity is evaluated independently with
 * bh_ffr_radiative_luminosity_code(), so hot low-rate flows retain the
 * quadratic radiatively inefficient branch. */
static double bh_ffr_cold_driving_luminosity_code(double mdot_h, double mdot_edd)
{
  if(mdot_h == 0)
    return 0.0;
  if(!(mdot_h > 0) || !(mdot_edd > 0))
    terminate("BH_FFR: invalid cold driving rates H=%g Edd=%g", mdot_h, mdot_edd);

  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  const double dotm_h = mdot_h / mdot_edd;

  if(dotm_h <= 1.0)
    return All.BHRadiativeEfficiency * mdot_h * c_internal * c_internal;

  const double l_edd = All.BHRadiativeEfficiency * mdot_edd * c_internal * c_internal;
  return l_edd * (1.0 + log(dotm_h));
}

static double bh_ffr_cold_wind_rate_code(double mdot_h, double mdot_edd)
{
  if(mdot_h == 0)
    return 0.0;

  const double l_drive_code = bh_ffr_cold_driving_luminosity_code(mdot_h, mdot_edd);
  const double l_drive_cgs = bh_ffr_power_code_to_cgs(l_drive_code);
  const double mdot_w_msun_yr =
      BH_FFR_GOFFORD_MDOT_NORM_MSUN_YR * pow(l_drive_cgs / BH_FFR_GOFFORD_L_NORM_CGS, BH_FFR_GOFFORD_MDOT_EXPONENT);
  const double mdot_w_cgs = mdot_w_msun_yr * SOLAR_MASS / SEC_PER_YEAR;
  const double mdot_w_code = bh_ffr_rate_cgs_to_code(mdot_w_cgs);

  if(!isfinite(mdot_w_code) || mdot_w_code < 0)
    terminate("BH_FFR: invalid Gofford cold-wind rate=%g", mdot_w_code);

  return mdot_w_code;
}

static double bh_ffr_cold_horizon_rate_code(double mdot_processed, double mdot_edd)
{
  if(!isfinite(mdot_processed) || mdot_processed < 0 || !isfinite(mdot_edd) || mdot_edd < 0)
    terminate("BH_FFR: invalid cold partition rates proc=%g Edd=%g", mdot_processed, mdot_edd);
  if(mdot_processed == 0)
    return 0.0;
  if(!(mdot_edd > 0))
    terminate("BH_FFR: positive cold processed rate=%g with non-positive Eddington rate", mdot_processed);

  /* F(H)=H+Mdot_w[L(H)]-Mdot_processed is monotonic, with F(0)<0 and
   * F(Mdot_processed)>=0.  Bisection is cheap, deterministic, and cannot
   * violate the mass bracket. */
  double lo = 0.0;
  double hi = mdot_processed;

  for(int iter = 0; iter < 96; iter++)
    {
      const double mid = 0.5 * (lo + hi);
      const double fmid = mid + bh_ffr_cold_wind_rate_code(mid, mdot_edd) - mdot_processed;

      if(fmid > 0)
        hi = mid;
      else
        lo = mid;

      if((hi - lo) <= 2.0e-13 * fmax(mdot_processed, DBL_MIN))
        break;
    }

  const double mdot_h = 0.5 * (lo + hi);
  if(!isfinite(mdot_h) || mdot_h < 0 || mdot_h > mdot_processed)
    terminate("BH_FFR: invalid cold root H=%g from processed=%g", mdot_h, mdot_processed);

  return mdot_h;
}

static void bh_ffr_hot_properties(double dotm_processed, double *retention, double *wind_velocity_code)
{
  if(!isfinite(dotm_processed) || dotm_processed < 0)
    terminate("BH_FFR: invalid hot processed Eddington ratio=%g", dotm_processed);
  if(!(All.BHRHotMaxInRs > 3.0))
    terminate("BH_FFR: BHRHotMaxInRs=%g must exceed 3", All.BHRHotMaxInRs);

  double rtr_rs = DBL_MAX;
  if(dotm_processed > 0)
    {
      const double q = BH_FFR_COLD_NOMINAL_EDD_RATIO / dotm_processed;
      if(q <= sqrt(DBL_MAX / 3.0))
        rtr_rs = 3.0 * q * q;
    }

  /* Hysteresis can briefly retain a hot label above the nominal cold
   * boundary, where the analytic Rtr would fall below 3 Rs.  Clamp the
   * inner radius at 3 Rs so f_H never exceeds unity. */
  const double r0_rs = fmax(3.0, fmin(rtr_rs, All.BHRHotMaxInRs));
  const double f_h = sqrt(3.0 / r0_rs);
  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  const double v_w = BH_FFR_HOT_WIND_VK_FACTOR * c_internal / sqrt(2.0 * r0_rs);

  if(!isfinite(f_h) || f_h < 0 || f_h > 1.0 || !isfinite(v_w) || v_w < 0)
    terminate("BH_FFR: invalid hot properties fH=%g vw=%g r0/Rs=%g", f_h, v_w, r0_rs);

  *retention = f_h;
  *wind_velocity_code = v_w;
}

static double bh_ffr_cold_wind_velocity_code(double mdot_h, double mdot_edd)
{
  if(mdot_h == 0)
    return 0.0;

  const double l_drive_cgs = bh_ffr_power_code_to_cgs(bh_ffr_cold_driving_luminosity_code(mdot_h, mdot_edd));
  double v_kms = BH_FFR_GOFFORD_V_NORM_KMS *
                 pow(l_drive_cgs / BH_FFR_GOFFORD_L_NORM_CGS, BH_FFR_GOFFORD_V_EXPONENT);
  v_kms = fmin(v_kms, BH_FFR_GOFFORD_V_MAX_KMS);

  const double v_code = v_kms * 1.0e5 / All.UnitVelocity_in_cm_per_s;
  if(!isfinite(v_code) || v_code < 0)
    terminate("BH_FFR: invalid cold wind velocity=%g", v_code);

  return v_code;
}

void bh_ffr_apply_inner_flow(int p, double processable_mass, double dt_code)
{
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid particle index=%d in inner-flow update", p);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for particle ID=%llu in inner-flow update", (unsigned long long)P[p].ID);

  if(!isfinite(processable_mass) || processable_mass < 0 || processable_mass > BHP[b].ReservoirMass ||
     !isfinite(dt_code) || dt_code < 0)
    terminate("BH_FFR: invalid inner-flow transaction dMproc=%g Md=%g dt=%g", processable_mass, BHP[b].ReservoirMass,
              dt_code);

  const double mdot_processed = BHP[b].MdotProcessed;
  const double mdot_edd = BHP[b].MdotEddington;

  /* An empty reservoir over a finite physical timestep is a valid no-op:
   * dt>0 with dMproc=0 must not terminate. The forbidden case is the
   * converse, a positive processed mass with zero elapsed time. */
  if(dt_code == 0 && processable_mass > 0)
    terminate("BH_FFR: positive inner-flow transaction dMproc=%g at zero dt", processable_mass);

  if(processable_mass == 0 && mdot_processed != 0)
    terminate("BH_FFR: zero processed mass with nonzero processed rate=%g at dt=%g", mdot_processed, dt_code);

  double mdot_h = 0.0;
  double wind_velocity_code = 0.0;

  if(mdot_processed > 0)
    {
      if(BHP[b].AccretionState == BH_FFR_STATE_COLD)
        {
          mdot_h = bh_ffr_cold_horizon_rate_code(mdot_processed, mdot_edd);
          wind_velocity_code = bh_ffr_cold_wind_velocity_code(mdot_h, mdot_edd);
        }
      else if(BHP[b].AccretionState == BH_FFR_STATE_ADIOS || BHP[b].AccretionState == BH_FFR_STATE_TRUNCATED)
        {
          double retention = 0.0;
          bh_ffr_hot_properties(BHP[b].ProcessedEddingtonRatio, &retention, &wind_velocity_code);
          mdot_h = retention * mdot_processed;
        }
      else
        terminate("BH_FFR: unphysical accretion state=%d in inner-flow update", BHP[b].AccretionState);
    }

  if(mdot_h < 0 && mdot_h > -1.0e-14 * fmax(1.0, mdot_processed))
    mdot_h = 0.0;
  if(mdot_h > mdot_processed && mdot_h < mdot_processed * (1.0 + 1.0e-13))
    mdot_h = mdot_processed;
  if(!isfinite(mdot_h) || mdot_h < 0 || mdot_h > mdot_processed)
    terminate("BH_FFR: invalid horizon rate H=%g processed=%g", mdot_h, mdot_processed);

  const double mdot_w = mdot_processed - mdot_h;
  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  const double lbol = bh_ffr_radiative_luminosity_code(mdot_h, mdot_edd);
  const double pwind = 0.5 * mdot_w * wind_velocity_code * wind_velocity_code;
  const double pjet =
      (BHP[b].AccretionState == BH_FFR_STATE_COLD) ? 0.0 : BH_FFR_HOT_JET_EFFICIENCY * mdot_h * c_internal * c_internal;

  double epsilon_rad = 0.0;
  if(mdot_h > 0)
    epsilon_rad = lbol / (mdot_h * c_internal * c_internal);

  if(!isfinite(mdot_w) || mdot_w < 0 || !isfinite(lbol) || lbol < 0 || !isfinite(pwind) || pwind < 0 ||
     !isfinite(pjet) || pjet < 0 || !isfinite(epsilon_rad) || epsilon_rad < 0 || epsilon_rad >= 1.0)
    terminate("BH_FFR: invalid Iteration-6 energetics for particle ID=%llu", (unsigned long long)P[p].ID);

  double delta_h = 0.0;
  if(mdot_processed > 0)
    delta_h = processable_mass * (mdot_h / mdot_processed);
  if(delta_h < 0 && delta_h > -1.0e-14 * fmax(1.0, processable_mass))
    delta_h = 0.0;
  if(delta_h > processable_mass && delta_h < processable_mass * (1.0 + 1.0e-13))
    delta_h = processable_mass;
  if(!isfinite(delta_h) || delta_h < 0 || delta_h > processable_mass)
    terminate("BH_FFR: invalid horizon mass dMH=%g from dMproc=%g", delta_h, processable_mass);

  const double delta_w = processable_mass - delta_h;
  const double delta_bh = (1.0 - epsilon_rad) * delta_h;
  const double reservoir_before = BHP[b].ReservoirMass;

  BHP[b].ReservoirMass = reservoir_before - processable_mass;
  if(BHP[b].ReservoirMass < 0 && BHP[b].ReservoirMass > -1.0e-13 * fmax(1.0, reservoir_before))
    BHP[b].ReservoirMass = 0.0;

  /* Removing a representative fraction of the unresolved reservoir should
   * preserve its coherence fraction rather than artificially increasing it. */
  const double coherence_scale = (reservoir_before > 0) ? BHP[b].ReservoirMass / reservoir_before : 0.0;
  for(int k = 0; k < 3; k++)
    BHP[b].Coherence[k] *= coherence_scale;

  BHP[b].BHMass += delta_bh;
  BHP[b].WindMassBuffer += delta_w;
  BHP[b].WindMomentumBuffer += delta_w * wind_velocity_code;
  BHP[b].WindEnergyBuffer += pwind * dt_code;
  BHP[b].JetEnergyBuffer += pjet * dt_code;

  BHP[b].MdotHorizon = mdot_h;
  BHP[b].MdotWind = mdot_w;
  BHP[b].BolometricLuminosity = lbol;
  BHP[b].WindPower = pwind;
  BHP[b].JetPower = pjet;

  /* Baseline mass-energy convention: BHMass retains horizon rest mass minus
   * radiated mass-energy. Mechanical powers are effective feedback energies
   * and are not subtracted a second time from BHMass. */
  P[p].Mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;

  const double partition_scale = fmax(1.0, processable_mass);
  if(fabs((delta_h + delta_w) - processable_mass) > 2.0e-12 * partition_scale)
    terminate("BH_FFR: Iteration-6 mass partition failed dMH+dMw=%g dMproc=%g", delta_h + delta_w, processable_mass);

  const double expected_dyn = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
  if(!isfinite(P[p].Mass) || P[p].Mass < 0 || fabs(P[p].Mass - expected_dyn) > 1.0e-12 * fmax(1.0, expected_dyn))
    terminate("BH_FFR: Iteration-6 dynamical-mass ledger failed for particle ID=%llu", (unsigned long long)P[p].ID);
}

void bh_ffr_inner_self_test(void)
{
  if(!(All.BHRHotMaxInRs > 3.0))
    terminate("BH_FFR: inner-flow self-test requires BHRHotMaxInRs>3");

  const double boundary_dotm = BH_FFR_COLD_NOMINAL_EDD_RATIO * sqrt(3.0 / All.BHRHotMaxInRs);
  double retention = 0.0, vhot = 0.0;
  bh_ffr_hot_properties(boundary_dotm, &retention, &vhot);
  const double expected_retention = sqrt(3.0 / All.BHRHotMaxInRs);
  if(fabs(retention - expected_retention) > 2.0e-13 * expected_retention)
    terminate("BH_FFR: inner-flow self-test failed hot retention got=%g expected=%g", retention, expected_retention);

  const double mbh_code = 1.0e5 * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;
  const double mdot_edd = bh_ffr_eddington_rate_code(mbh_code);
  if(!(mdot_edd > 0))
    terminate("BH_FFR: inner-flow self-test invalid Eddington rate");

  const double cold_proc = 10.0 * mdot_edd;
  const double cold_h = bh_ffr_cold_horizon_rate_code(cold_proc, mdot_edd);
  const double cold_w = bh_ffr_cold_wind_rate_code(cold_h, mdot_edd);
  if(fabs((cold_h + cold_w) - cold_proc) > 2.0e-11 * cold_proc)
    terminate("BH_FFR: inner-flow self-test failed cold mass closure H=%g W=%g proc=%g", cold_h, cold_w, cold_proc);
  if(!(cold_h > mdot_edd))
    terminate("BH_FFR: inner-flow self-test accidentally imposed an Eddington mass cap H/Edd=%g", cold_h / mdot_edd);

  const double l1_low = bh_ffr_radiative_luminosity_code(mdot_edd, mdot_edd);
  const double c_internal = CLIGHT / All.UnitVelocity_in_cm_per_s;
  const double l_edd_code = All.BHRadiativeEfficiency * mdot_edd * c_internal * c_internal;
  if(fabs(l1_low - l_edd_code) > 2.0e-13 * l_edd_code)
    terminate("BH_FFR: inner-flow self-test failed luminosity continuity at dotm_H=1");

  const double l_low_a = bh_ffr_radiative_luminosity_code(0.01 * mdot_edd, mdot_edd);
  const double l_low_b = bh_ffr_radiative_luminosity_code(0.02 * mdot_edd, mdot_edd);
  if(fabs(l_low_b / l_low_a - 4.0) > 2.0e-12)
    terminate("BH_FFR: inner-flow self-test failed low-rate quadratic luminosity scaling");

  const double pjet = BH_FFR_HOT_JET_EFFICIENCY * mdot_edd * c_internal * c_internal;
  if(!(vhot > 0) || !(pjet > 0))
    terminate("BH_FFR: inner-flow self-test failed hot wind/jet energetics");
}
