#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

static int bh_ffr_compact_index_from_particle(int p, const char *where)
{
  if(p < 0 || p >= NumPart)
    terminate("BH_FFR: particle index %d outside NumPart=%d in %s", p, NumPart, where);
  if(P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: non-Type-%d particle ID=%llu passed to %s", BH_FFR_PARTICLE_TYPE,
              (unsigned long long)P[p].ID, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  return b;
}

/* The PDF requires one beginning-of-transaction reservoir time to control
 * drainage, jet-axis relaxation, and the next reservoir-accuracy timestep.
 * Keep this last-transaction value in a small transient ID-addressed cache:
 * it is purely timestep bookkeeping and does not need another persistent BHP
 * field even though the benchmark BHP layout is explicitly versioned.
 *
 * The post-transaction BH/reservoir masses are a signature: a subsequent BH
 * merger invalidates the cached value automatically, while wind/jet packet
 * release does not.  Stale slots are recycled on the next synchronization. */
struct bh_ffr_frozen_disk_time_entry
{
  MyIDType ID;
  integertime Ti;
  double DiskTimeMyr;
  double PostBHMass;
  double PostReservoirMass;
};

static struct bh_ffr_frozen_disk_time_entry *FrozenDiskTimeCache;
static int FrozenDiskTimeCount;
static int FrozenDiskTimeCapacity;

static void bh_ffr_store_frozen_disk_time(int p, int b, double disk_time_myr)
{
  if(!isfinite(disk_time_myr) || !(disk_time_myr > 0))
    terminate("BH_FFR: invalid frozen reservoir time=%g Myr for ID=%llu", disk_time_myr,
              (unsigned long long)P[p].ID);

  int slot = -1;
  int stale = -1;
  for(int i = 0; i < FrozenDiskTimeCount; i++)
    {
      if(FrozenDiskTimeCache[i].ID == P[p].ID)
        {
          slot = i;
          break;
        }
      if(stale < 0 && FrozenDiskTimeCache[i].Ti != All.Ti_Current)
        stale = i;
    }

  if(slot < 0)
    slot = stale;

  if(slot < 0)
    {
      if(FrozenDiskTimeCount == FrozenDiskTimeCapacity)
        {
          const int new_capacity = FrozenDiskTimeCapacity > 0 ? 2 * FrozenDiskTimeCapacity : 16;
          void *tmp = realloc(FrozenDiskTimeCache, (size_t)new_capacity * sizeof(*FrozenDiskTimeCache));
          if(tmp == NULL)
            terminate("BH_FFR: failed to grow frozen-reservoir-time cache to %d entries", new_capacity);
          FrozenDiskTimeCache = (struct bh_ffr_frozen_disk_time_entry *)tmp;
          FrozenDiskTimeCapacity = new_capacity;
        }
      slot = FrozenDiskTimeCount++;
    }

  FrozenDiskTimeCache[slot].ID = P[p].ID;
  FrozenDiskTimeCache[slot].Ti = All.Ti_Current;
  FrozenDiskTimeCache[slot].DiskTimeMyr = disk_time_myr;
  FrozenDiskTimeCache[slot].PostBHMass = BHP[b].BHMass;
  FrozenDiskTimeCache[slot].PostReservoirMass = BHP[b].ReservoirMass;
}

static double bh_ffr_timestep_reservoir_time_myr(int p, int b, int *used_frozen)
{
  for(int i = 0; i < FrozenDiskTimeCount; i++)
    if(FrozenDiskTimeCache[i].ID == P[p].ID && FrozenDiskTimeCache[i].Ti == All.Ti_Current &&
       FrozenDiskTimeCache[i].PostBHMass == BHP[b].BHMass &&
       FrozenDiskTimeCache[i].PostReservoirMass == BHP[b].ReservoirMass)
      {
        if(used_frozen != NULL)
          *used_frozen = 1;
        return FrozenDiskTimeCache[i].DiskTimeMyr;
      }

  /* Startup/newly seeded BHs and post-merger survivors can legitimately have
   * no matching transaction cache.  In those cases limit from the current
   * state; ordinary processed BHs must take the frozen branch above. */
  if(used_frozen != NULL)
    *used_frozen = 0;
  return bh_ffr_reservoir_timescale_myr(BHP[b].BHMass, BHP[b].ReservoirMass);
}

double bh_ffr_integer_interval_to_physical_myr(integertime ti0, integertime ti1)
{
  if(ti0 < 0 || ti1 < ti0 || ti1 > TIMEBASE)
    terminate("BH_FFR: invalid integer-time interval [%lld,%lld] with TIMEBASE=%lld", (long long)ti0, (long long)ti1,
              (long long)TIMEBASE);

  if(ti0 == ti1)
    return 0.0;

  double x0, x1;
  if(All.ComovingIntegrationOn)
    {
      x0 = All.TimeBegin * exp((double)ti0 * All.Timebase_interval);
      x1 = All.TimeBegin * exp((double)ti1 * All.Timebase_interval);
    }
  else
    {
      x0 = All.TimeBegin + (double)ti0 * All.Timebase_interval;
      x1 = All.TimeBegin + (double)ti1 * All.Timebase_interval;
    }

  const double dt_myr = 1000.0 * get_time_difference_in_Gyr(x0, x1);
  if(!isfinite(dt_myr) || dt_myr < 0)
    terminate("BH_FFR: invalid physical timestep %g Myr from integer interval [%lld,%lld]", dt_myr, (long long)ti0,
              (long long)ti1);

  return dt_myr;
}

double bh_ffr_integer_interval_to_physical_code_time(integertime ti0, integertime ti1)
{
  const double dt_myr = bh_ffr_integer_interval_to_physical_myr(ti0, ti1);
  if(dt_myr == 0)
    return 0.0;

  if(!(All.UnitTime_in_s > 0) || !(All.HubbleParam > 0))
    terminate("BH_FFR: invalid unit conversion UnitTime=%g HubbleParam=%g", All.UnitTime_in_s, All.HubbleParam);

  const double dt = dt_myr * SEC_PER_MEGAYEAR * All.HubbleParam / All.UnitTime_in_s;
  if(!isfinite(dt) || dt < 0)
    terminate("BH_FFR: invalid physical code-time interval %g", dt);

  return dt;
}

double bh_ffr_get_elapsed_time_myr(int p)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_get_elapsed_time_myr");

  if(BHP[b].LastProcessedTi > All.Ti_Current)
    terminate("BH_FFR: LastProcessedTi=%lld exceeds Ti_Current=%lld for particle ID=%llu", (long long)BHP[b].LastProcessedTi,
              (long long)All.Ti_Current, (unsigned long long)P[p].ID);

  return bh_ffr_integer_interval_to_physical_myr(BHP[b].LastProcessedTi, All.Ti_Current);
}

double bh_ffr_get_elapsed_time_code_time(int p)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_get_elapsed_time_code_time");

  if(BHP[b].LastProcessedTi > All.Ti_Current)
    terminate("BH_FFR: LastProcessedTi=%lld exceeds Ti_Current=%lld for particle ID=%llu", (long long)BHP[b].LastProcessedTi,
              (long long)All.Ti_Current, (unsigned long long)P[p].ID);

  return bh_ffr_integer_interval_to_physical_code_time(BHP[b].LastProcessedTi, All.Ti_Current);
}

void bh_ffr_set_min_neighbour_timebin(int p, int timebin)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_set_min_neighbour_timebin");

  if(timebin < -1 || timebin == 0 || timebin >= TIMEBINS)
    terminate("BH_FFR: invalid neighbour hydro timebin=%d for particle ID=%llu", timebin, (unsigned long long)P[p].ID);

  BHP[b].MinNeighbourHydroTimeBin = timebin;
}

static double bh_ffr_code_time_to_myr(double dt_code)
{
  if(!isfinite(dt_code) || dt_code < 0 || !(All.UnitTime_in_s > 0) || !(All.HubbleParam > 0))
    terminate("BH_FFR: invalid code-time to Myr conversion dt=%g", dt_code);

  return dt_code * All.UnitTime_in_s / (All.HubbleParam * SEC_PER_MEGAYEAR);
}

static integertime bh_ffr_limit_integer_step_by_physical_myr(integertime ti_step, double limit_myr)
{
  if(!isfinite(limit_myr) || !(limit_myr > 0) || ti_step <= 2)
    return ti_step;

  integertime limited = ti_step;

  while(limited > 2)
    {
      integertime ti1 = All.Ti_Current + limited;
      if(ti1 > TIMEBASE)
        ti1 = TIMEBASE;

      const double dt_myr = bh_ffr_integer_interval_to_physical_myr(All.Ti_Current, ti1);
      if(dt_myr <= limit_myr)
        break;

      limited >>= 1;
      if(limited < 2)
        limited = 2;
    }

  return limited;
}

static double bh_ffr_feedback_channel_limit_myr(double threshold,
                                                double power)
{
  if(!(threshold > 0) || !(power > 0))
    return HUGE_VAL;

  const double tth_myr = bh_ffr_code_time_to_myr(threshold / power);
  const double limit =
      BH_FFR_FEEDBACK_TARGET_QUANTA_PER_STEP * tth_myr;

  if(!isfinite(limit) || !(limit > 0))
    terminate("BH_FFR: invalid feedback cadence limit Eth=%g P=%g limit=%g Myr",
              threshold, power, limit);
  return limit;
}

static int bh_ffr_feedback_target_hydro_timebin_core(int p,
                                                     double threshold,
                                                     double power)
{
  if(!(threshold > 0) || !(power > 0))
    return -1;

  const int grav_bin = P[p].TimeBinGrav;
  if(grav_bin <= 0 || grav_bin >= TIMEBINS)
    terminate("BH_FFR: invalid BH gravity timebin=%d for feedback wake ID=%llu",
              grav_bin, (unsigned long long)P[p].ID);

  integertime ti_step = ((integertime)1) << grav_bin;
  const double limit_myr =
      bh_ffr_feedback_channel_limit_myr(threshold, power);
  ti_step =
      bh_ffr_limit_integer_step_by_physical_myr(ti_step, limit_myr);

  const int bin = get_timestep_bin(ti_step);
  if(bin <= 0 || bin >= TIMEBINS)
    terminate("BH_FFR: invalid requested feedback hydro timebin=%d ti=%lld ID=%llu",
              bin, (long long)ti_step, (unsigned long long)P[p].ID);
  return bin;
}

int bh_ffr_feedback_wind_target_hydro_timebin(int p)
{
  const int b =
      bh_ffr_compact_index_from_particle(
          p, "bh_ffr_feedback_wind_target_hydro_timebin");
  return bh_ffr_feedback_target_hydro_timebin_core(
      p, BHP[b].WindThresholdEnergy, BHP[b].WindPower);
}

int bh_ffr_feedback_jet_target_hydro_timebin(int p)
{
  const int b =
      bh_ffr_compact_index_from_particle(
          p, "bh_ffr_feedback_jet_target_hydro_timebin");
  return bh_ffr_feedback_target_hydro_timebin_core(
      p, BHP[b].JetThresholdEnergy, BHP[b].JetPower);
}

/* Shared guard-zone cadence for gas inside the maximum MACER receiver sphere.
 * Start from the actual BH gravity timebin and only shorten it further if the
 * wind or jet threshold/power cadence requires a smaller synchronized step.
 * Gas already on a finer hydro step is never lengthened by the wake manager. */
int bh_ffr_feedback_sync_target_hydro_timebin(int p)
{
  const int b =
      bh_ffr_compact_index_from_particle(
          p, "bh_ffr_feedback_sync_target_hydro_timebin");

  const int grav_bin = P[p].TimeBinGrav;
  if(grav_bin <= 0 || grav_bin >= TIMEBINS)
    terminate("BH_FFR: invalid BH gravity timebin=%d for feedback sync guard ID=%llu",
              grav_bin, (unsigned long long)P[p].ID);

  integertime ti_step = ((integertime)1) << grav_bin;

  const double wind_limit_myr =
      bh_ffr_feedback_channel_limit_myr(
          BHP[b].WindThresholdEnergy, BHP[b].WindPower);
  const double jet_limit_myr =
      bh_ffr_feedback_channel_limit_myr(
          BHP[b].JetThresholdEnergy, BHP[b].JetPower);

  ti_step =
      bh_ffr_limit_integer_step_by_physical_myr(ti_step, wind_limit_myr);
  ti_step =
      bh_ffr_limit_integer_step_by_physical_myr(ti_step, jet_limit_myr);

  const int bin = get_timestep_bin(ti_step);
  if(bin <= 0 || bin >= TIMEBINS)
    terminate("BH_FFR: invalid feedback sync-guard timebin=%d ti=%lld ID=%llu",
              bin, (long long)ti_step, (unsigned long long)P[p].ID);

  return bin;
}

static int *BHFFRFeedbackWakeBin = NULL;
static unsigned char *BHFFRFeedbackHydroActiveNow = NULL;

void bh_ffr_feedback_wakeup_begin(void)
{
  if(BHFFRFeedbackWakeBin != NULL || BHFFRFeedbackHydroActiveNow != NULL)
    terminate("BH_FFR: feedback wakeup transaction already open");

  const size_t ngas_alloc = (size_t)(NumGas > 0 ? NumGas : 1);

  BHFFRFeedbackWakeBin =
      (int *)malloc(ngas_alloc * sizeof(int));
  BHFFRFeedbackHydroActiveNow =
      (unsigned char *)malloc(ngas_alloc * sizeof(unsigned char));

  if(BHFFRFeedbackWakeBin == NULL || BHFFRFeedbackHydroActiveNow == NULL)
    terminate("BH_FFR: failed to allocate feedback wakeup/activity state");

  for(int gas = 0; gas < NumGas; gas++)
    {
      BHFFRFeedbackWakeBin[gas] = TIMEBINS;
      BHFFRFeedbackHydroActiveNow[gas] = 0;
    }

  /* Snapshot the actual current-step hydro active list before feedback.  The
   * active cells may already have been assigned a different TimeBinHydro for
   * their NEXT step, so TimeBinSynchronized[P[i].TimeBinHydro] is not a
   * reliable indicator of whether they just completed hydro at Ti_Current. */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      const int gas = TimeBinsHydro.ActiveParticleList[idx];
      if(gas < 0)
        continue;
      if(gas >= NumGas || P[gas].Type != 0)
        terminate("BH_FFR: invalid hydro-active receiver index=%d NumGas=%d",
                  gas, NumGas);
      if(P[gas].ID == 0 || P[gas].Mass == 0)
        continue;
      BHFFRFeedbackHydroActiveNow[gas] = 1;
    }
}

int bh_ffr_feedback_gas_is_hydro_active_now(int gas_index)
{
  if(BHFFRFeedbackHydroActiveNow == NULL)
    terminate("BH_FFR: feedback hydro-active query outside feedback transaction");
  if(gas_index < 0 || gas_index >= NumGas)
    terminate("BH_FFR: hydro-active gas index=%d outside NumGas=%d",
              gas_index, NumGas);

  return BHFFRFeedbackHydroActiveNow[gas_index] != 0;
}

void bh_ffr_feedback_request_wakeup(int gas_index, int target_timebin)
{
  if(BHFFRFeedbackWakeBin == NULL || target_timebin < 0)
    return;
  if(gas_index < 0 || gas_index >= NumGas)
    terminate("BH_FFR: wakeup gas index=%d outside NumGas=%d",
              gas_index, NumGas);
  if(target_timebin <= 0 || target_timebin >= TIMEBINS)
    terminate("BH_FFR: invalid feedback wakeup target bin=%d gas=%d",
              target_timebin, gas_index);

  int old = BHFFRFeedbackWakeBin[gas_index];
  while(target_timebin < old)
    {
      const int seen =
          __sync_val_compare_and_swap(&BHFFRFeedbackWakeBin[gas_index],
                                      old, target_timebin);
      if(seen == old)
        break;
      old = seen;
    }
}

void bh_ffr_feedback_wakeup_apply(void)
{
  if(BHFFRFeedbackWakeBin == NULL)
    return;

  int requested = 0;
  int moved = 0;
  int min_new_bin = TIMEBINS;
  int max_old_bin = -1;

  for(int gas = 0; gas < NumGas; gas++)
    {
      int target = BHFFRFeedbackWakeBin[gas];
      if(target >= TIMEBINS)
        continue;
      requested++;

      const int encoded_old = P[gas].TimeBinHydro;
      int old = encoded_old;
      if(old < 0)
        old = -old - 1;

      if(old <= 0 || old >= TIMEBINS)
        terminate("BH_FFR: invalid gas hydro timebin=%d decoded=%d gasID=%llu during wakeup",
                  encoded_old, old, (unsigned long long)P[gas].ID);

      while(target > 1 && !TimeBinSynchronized[target])
        target--;

      if(target <= 0 || target >= TIMEBINS)
        terminate("BH_FFR: failed to find synchronized feedback wakeup bin for gasID=%llu",
                  (unsigned long long)P[gas].ID);

      if(target >= old)
        continue;

      timebin_move_particle(&TimeBinsHydro, gas, old, target);
      P[gas].TimeBinHydro =
          encoded_old < 0 ? -target - 1 : target;

      moved++;
      if(target < min_new_bin)
        min_new_bin = target;
      if(old > max_old_bin)
        max_old_bin = old;
    }

  printf("BH_FFR: feedback wakeup task=%d requested=%d moved=%d "
         "minNewBin=%d maxOldBin=%d\n",
         ThisTask, requested, moved,
         moved > 0 ? min_new_bin : -1,
         moved > 0 ? max_old_bin : -1);
  fflush(stdout);

  free(BHFFRFeedbackWakeBin);
  free(BHFFRFeedbackHydroActiveNow);
  BHFFRFeedbackWakeBin = NULL;
  BHFFRFeedbackHydroActiveNow = NULL;
}

integertime bh_ffr_limit_gravity_timestep(int p, integertime ti_step)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_limit_gravity_timestep");
  integertime limited = ti_step;
  const integertime raw_step = ti_step;

  /* Keep the BH no coarser than the fastest gas inside the common accretion
   * aperture.  Feedback itself is coupled only when receiving gas is hydro
   * active, so feedback-buffer thresholds must not force still finer BH-only
   * subcycles. */
  const int gas_bin = BHP[b].MinNeighbourHydroTimeBin;
  if(gas_bin >= 0)
    {
      if(gas_bin == 0 || gas_bin >= TIMEBINS)
        terminate("BH_FFR: corrupt neighbour hydro timebin=%d for particle ID=%llu",
                  gas_bin, (unsigned long long)P[p].ID);

      const integertime gas_ti_step = ((integertime)1) << gas_bin;
      if(gas_ti_step < limited)
        limited = gas_ti_step;
    }

  double dt_limit_myr = HUGE_VAL;
  double reservoir_limit_myr = HUGE_VAL;
  double wind_limit_myr = HUGE_VAL;
  double jet_limit_myr = HUGE_VAL;
  double tng_limit_myr = HUGE_VAL;
  int used_frozen_disk_time = 0;

  /* Reservoir evolution and Chandrasekhar drag are genuine internal dynamical
   * timescales and may require the BH to subcycle independently of hydro. */
  if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_RESERVOIR)
    {
      const double disk_time_myr =
          bh_ffr_timestep_reservoir_time_myr(p, b, &used_frozen_disk_time);
      if(disk_time_myr < DBL_MAX / All.BHInternalTimestepFactor)
        reservoir_limit_myr = All.BHInternalTimestepFactor * disk_time_myr;
      if(reservoir_limit_myr < dt_limit_myr)
        dt_limit_myr = reservoir_limit_myr;

      const double df_tcode =
          bh_ffr_get_cached_dynamical_friction_timescale_code(p);
      if(isfinite(df_tcode) && df_tcode > 0)
        {
          const double df_limit_myr =
              All.BHInternalTimestepFactor * bh_ffr_code_time_to_myr(df_tcode);
          if(!isfinite(df_limit_myr) || !(df_limit_myr > 0))
            terminate("BH_FFR: invalid dynamical-friction timestep limit=%g Myr for ID=%llu",
                      df_limit_myr, (unsigned long long)P[p].ID);
          if(df_limit_myr < dt_limit_myr)
            dt_limit_myr = df_limit_myr;
        }
    }

  /* Use one shared physical feedback cadence for source generation and
   * receiver response.  The BH source step and the deferred receiver-gas wake
   * both use Eth/P, preventing either side from outrunning the other. */
  if(All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_MACER)
    {
      wind_limit_myr =
          bh_ffr_feedback_channel_limit_myr(
              BHP[b].WindThresholdEnergy, BHP[b].WindPower);
      jet_limit_myr =
          bh_ffr_feedback_channel_limit_myr(
              BHP[b].JetThresholdEnergy, BHP[b].JetPower);

      if(wind_limit_myr < dt_limit_myr)
        dt_limit_myr = wind_limit_myr;
      if(jet_limit_myr < dt_limit_myr)
        dt_limit_myr = jet_limit_myr;
    }
  else if(All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_TNG &&
          BHP[b].TNGFeedbackMode == BH_BENCHMARK_TNG_MODE_KINETIC &&
          BHP[b].TNGKineticThresholdEnergy > 0 && BHP[b].TNGFeedbackPower > 0)
    {
      const double tcode =
          BHP[b].TNGKineticThresholdEnergy / BHP[b].TNGFeedbackPower;
      tng_limit_myr =
          All.BHInternalTimestepFactor * bh_ffr_code_time_to_myr(tcode);
      if(!isfinite(tng_limit_myr) || !(tng_limit_myr > 0))
        terminate("BH_TNG: invalid kinetic threshold time=%g Myr for ID=%llu",
                  tng_limit_myr, (unsigned long long)P[p].ID);
    }

  limited = bh_ffr_limit_integer_step_by_physical_myr(limited, dt_limit_myr);
  const integertime accuracy_limited = limited;

  /* Backlog remains a diagnostic only.  Buffered feedback is discharged when
   * synchronized hydro targets satisfy the coupling criteria; making the BH
   * timestep even shorter cannot create such targets. */
  int backlog = 0;
  if(All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_MACER)
    {
      if(BHP[b].WindThresholdEnergy > 0 &&
         BHP[b].WindEnergyBuffer >= BHP[b].WindThresholdEnergy)
        backlog = 1;
      if(BHP[b].JetThresholdEnergy > 0 &&
         BHP[b].JetEnergyBuffer >= BHP[b].JetThresholdEnergy)
        backlog = 1;
    }
  else if(All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_TNG)
    {
      if(BHP[b].TNGKineticThresholdEnergy > 0 &&
         BHP[b].TNGKineticEnergyBuffer >= BHP[b].TNGKineticThresholdEnergy)
        backlog = 1;
      if(BHP[b].TNGThermalEnergyBuffer > 0)
        backlog = 1;
    }

  const integertime backlog_limit = limited;
  const int backlog_applied = 0;

  if(limited < raw_step)
    {
      const double wind_ratio =
          BHP[b].WindThresholdEnergy > 0
              ? BHP[b].WindEnergyBuffer / BHP[b].WindThresholdEnergy
              : 0.0;
      const double jet_ratio =
          BHP[b].JetThresholdEnergy > 0
              ? BHP[b].JetEnergyBuffer / BHP[b].JetThresholdEnergy
              : 0.0;
      const double tng_ratio =
          BHP[b].TNGKineticThresholdEnergy > 0
              ? BHP[b].TNGKineticEnergyBuffer /
                    BHP[b].TNGKineticThresholdEnergy
              : 0.0;

      printf("BH_FFR: timestep limit ID=%llu task=%d feedback=%d raw=%lld limited=%lld "
             "gasbin=%d fint=%g dtresMyr=%g dtwindMyr=%g dtjetMyr=%g dtTNGMyr=%g "
             "accuracyLimited=%lld backlogLimit=%lld windBacklog=%g jetBacklog=%g "
             "tngBacklog=%g backlog=%d backlogApplied=%d tdFrozen=%d\n",
             (unsigned long long)P[p].ID, ThisTask,
             All.BHBenchmarkFeedbackModel, (long long)raw_step,
             (long long)limited, gas_bin, All.BHInternalTimestepFactor,
             reservoir_limit_myr, wind_limit_myr, jet_limit_myr,
             tng_limit_myr, (long long)accuracy_limited,
             (long long)backlog_limit, wind_ratio, jet_ratio, tng_ratio,
             backlog, backlog_applied, used_frozen_disk_time);
      fflush(stdout);
    }

  return limited;
}

/* One-Myr physical seed-quiescence experiment (this branch only).
 * ColdBlendWeight is the unused restart-persistent Iteration-5 slot. We
 * temporarily store a_seed here without changing the native restart ABI.
 * In the pristine post-seed z=20.8837 restart this slot is zero, so its
 * first BH step establishes the epoch. Newly seeded BHs do likewise. */
#define BH_FFR_SEED_QUIET_MYR 1.0
int bh_ffr_seed_quiescent(int p)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_seed_quiescent");

  if(!All.ComovingIntegrationOn || !(All.Time > 0.0) || !(All.Time < 1.0))
    terminate("BH_FFR: 1-Myr seed delay requires cosmological 0<a<1 (a=%g)", All.Time);

  if(BHP[b].ColdBlendWeight == 0.0)
    {
      BHP[b].ColdBlendWeight = All.Time;
      printf("BH_FFR: seed delay started ID=%llu a_seed=%.17g durationMyr=%g\n",
             (unsigned long long)P[p].ID, All.Time, BH_FFR_SEED_QUIET_MYR);
      fflush(stdout);
    }

  const double a_seed = BHP[b].ColdBlendWeight;
  if(!isfinite(a_seed) || !(a_seed > 0.0) || a_seed > All.Time)
    terminate("BH_FFR: invalid seed epoch ID=%llu a_seed=%g a=%g",
              (unsigned long long)P[p].ID, a_seed, All.Time);

  const double age_myr = 1000.0 * get_time_difference_in_Gyr(a_seed, All.Time);
  if(!isfinite(age_myr) || age_myr < 0.0)
    terminate("BH_FFR: invalid seed age ID=%llu ageMyr=%g",
              (unsigned long long)P[p].ID, age_myr);

  return age_myr < BH_FFR_SEED_QUIET_MYR;
}

void bh_ffr_step(void)
{
  /* Particle array indices can change after domain/FoF reordering. */
  bh_ffr_build_active_list();

  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

  static int self_tests_done = 0;
  if(!self_tests_done)
    {
      bh_ffr_adaptive_radius_self_test();
      bh_ffr_reservoir_self_test();
      bh_ffr_inner_self_test();
      bh_ffr_feedback_self_test();
      bh_ffr_jet_self_test();
      bh_benchmark_tng_feedback_self_test();
      self_tests_done = 1;
    }

  /* The selected resolved estimator removes gas exactly once and stores both
   * its uncapped algebraic rate and the realized conservative capture rate. */
  bh_ffr_capture_resolved_gas();

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_step");

      if(P[p].Ti_Current != All.Ti_Current)
        terminate("BH_FFR: active particle ID=%llu is not drifted to Ti_Current=%lld",
                  (unsigned long long)P[p].ID, (long long)All.Ti_Current);

      const double dt_myr = bh_ffr_get_elapsed_time_myr(p);
      const double dt_code = bh_ffr_get_elapsed_time_code_time(p);
      if(BHP[b].LastProcessedTi != All.Ti_Current &&
         (!(dt_myr > 0) || !(dt_code > 0)))
        terminate("BH_FFR: non-positive elapsed physical timestep for active particle ID=%llu",
                  (unsigned long long)P[p].ID);

      const double operational_mass = BHP[b].BenchmarkMdotOperational * dt_code;
      if(!isfinite(operational_mass) || operational_mass < 0)
        terminate("BH_BENCHMARK: invalid operational mass increment=%g for ID=%llu",
                  operational_mass, (unsigned long long)P[p].ID);
      BHP[b].BenchmarkCumulativeOperationalMass += operational_mass;

      if(All.BHBenchmarkAccretionTarget == BH_BENCHMARK_TARGET_RESERVOIR)
        {
          const double bh_mass_before_inner = BHP[b].BHMass;

          /* Preserve the validated FFR-MACER orbital and reservoir transaction
           * exactly for every reservoir-backed benchmark. */
          bh_ffr_apply_cached_dynamical_friction(p, dt_code);

          const double disk_time_myr =
              bh_ffr_reservoir_timescale_myr(BHP[b].BHMass,
                                              BHP[b].ReservoirMass);

          /* During seed quiescence neither horizon growth nor wind/jet
           * source generation is permitted; empty buffers remain empty. */
          if(!bh_ffr_seed_quiescent(p))
            bh_ffr_update_reservoir_state(p, dt_myr, dt_code, disk_time_myr);

          const double bh_growth = BHP[b].BHMass - bh_mass_before_inner;
          if(!isfinite(bh_growth) ||
             bh_growth < -2.0e-12 * fmax(fabs(bh_mass_before_inner), 1.0e-30))
            terminate("BH_BENCHMARK: invalid reservoir BH-mass growth=%g for ID=%llu",
                      bh_growth, (unsigned long long)P[p].ID);
          if(bh_growth > 0)
            BHP[b].BenchmarkCumulativeBHMassGrowth += bh_growth;

          if(!bh_ffr_seed_quiescent(p))
            {
              bh_ffr_update_jet_direction(p, dt_myr, disk_time_myr);
              bh_ffr_store_frozen_disk_time(p, b, disk_time_myr);
            }
        }
      else
        {
          /* Direct benchmark capture bypasses the FFR reservoir and MACER
           * inner-flow bookkeeping. MACER buffers must therefore remain empty;
           * the separate TNG buffers are allowed when TNG feedback is selected. */
          if(BHP[b].ReservoirMass != 0 || BHP[b].WindMassBuffer != 0 ||
             BHP[b].WindMomentumBuffer != 0 ||
             BHP[b].WindEnergyBuffer != 0 || BHP[b].JetEnergyBuffer != 0)
            terminate("BH_BENCHMARK: direct backend found non-empty MACER state for ID=%llu",
                      (unsigned long long)P[p].ID);
        }

      printf("BH_BENCHMARK: mass-ledger ID=%llu task=%d dtCode=%.17g "
             "dMop=%.17g cumOperational=%.17g cumRealized=%.17g cumBHGrowth=%.17g\n",
             (unsigned long long)P[p].ID, ThisTask, dt_code, operational_mass,
             BHP[b].BenchmarkCumulativeOperationalMass,
             BHP[b].BenchmarkCumulativeRealizedMass,
             BHP[b].BenchmarkCumulativeBHMassGrowth);
      fflush(stdout);
    }

  /* TNG feedback energy is generated from this transaction before
   * LastProcessedTi is advanced. Its mode is selected from the uncapped
   * estimator/Eddington ratio, while energy generation uses the selected
   * Eddington-limited operational accretion rate as in the TNG source model.
   * The realized finite-step gas sink remains a separate conservation
   * diagnostic and converges to that operational rate with timestep. */
  if(All.BHBenchmarkFeedbackModel == BH_BENCHMARK_FEEDBACK_TNG)
    bh_benchmark_tng_feedback_accumulate();

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_step-commit-time");
      BHP[b].LastProcessedTi = All.Ti_Current;
    }

  /* Merge after source accumulation but before injection. Persistent feedback
   * reservoirs are combined conservatively by the merger routine, so a pair
   * cannot emit two spatially overlapping events immediately before coalescing. */
  bh_ffr_merge_close_black_holes();

  /* No packet injection while every BH remains inside its seed quiet phase.
   * The experiment uses one BH; this also avoids discharging pre-existing
   * buffers from its pristine checkpoint before the delay expires. */
  int local_ready_bhs = 0, global_ready_bhs = 0;
  for(int n = 0; n < NumActiveBHFFR; n++)
    if(!bh_ffr_seed_quiescent(BHFFRActiveParticleList[n]))
      local_ready_bhs++;
  MPI_Allreduce(&local_ready_bhs, &global_ready_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_ready_bhs == 0)
    return;

  switch(All.BHBenchmarkFeedbackModel)
    {
      case BH_BENCHMARK_FEEDBACK_NONE:
        break;

      case BH_BENCHMARK_FEEDBACK_TNG:
        bh_benchmark_tng_feedback_inject();
        break;

      case BH_BENCHMARK_FEEDBACK_MACER:
        bh_ffr_feedback_wakeup_begin();
        bh_ffr_inject_wind_feedback();
        bh_ffr_inject_jet_feedback();
        bh_ffr_feedback_wakeup_apply();
        break;

      default:
        terminate("BH_BENCHMARK: unknown feedback model=%d",
                  All.BHBenchmarkFeedbackModel);
    }
}