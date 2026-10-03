#include <math.h>
#include <mpi.h>

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

integertime bh_ffr_limit_gravity_timestep(int p, integertime ti_step)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_limit_gravity_timestep");
  integertime limited = ti_step;
  const integertime raw_step = ti_step;

  /* Existing local-hydro synchronization constraint. */
  const int gas_bin = BHP[b].MinNeighbourHydroTimeBin;
  if(gas_bin >= 0)
    {
      if(gas_bin == 0 || gas_bin >= TIMEBINS)
        terminate("BH_FFR: corrupt neighbour hydro timebin=%d for particle ID=%llu", gas_bin, (unsigned long long)P[p].ID);

      const integertime gas_ti_step = ((integertime)1) << gas_bin;
      if(gas_ti_step < limited)
        limited = gas_ti_step;
    }

  /* Iteration 9 internal accuracy limits. The analytic reservoir/direction
   * maps are stable without this restriction; this controls coefficient and
   * burst-threshold evolution over one BH step.  Energy/power ratios are in
   * physical code time and are converted to Myr before timeline limiting. */
  const double reservoir_limit_myr = All.BHInternalTimestepFactor * All.BHDiskTimeMyr;
  double dt_limit_myr = reservoir_limit_myr;
  double wind_limit_myr = HUGE_VAL;
  double jet_limit_myr = HUGE_VAL;

  if(BHP[b].WindThresholdEnergy > 0 && BHP[b].WindPower > 0)
    {
      const double tcode = BHP[b].WindThresholdEnergy / BHP[b].WindPower;
      wind_limit_myr = All.BHInternalTimestepFactor * bh_ffr_code_time_to_myr(tcode);
      if(!isfinite(wind_limit_myr) || !(wind_limit_myr > 0))
        terminate("BH_FFR: invalid wind timestep limit=%g Myr for ID=%llu", wind_limit_myr,
                  (unsigned long long)P[p].ID);
      if(wind_limit_myr < dt_limit_myr)
        dt_limit_myr = wind_limit_myr;
    }

  if(BHP[b].JetThresholdEnergy > 0 && BHP[b].JetPower > 0)
    {
      const double tcode = BHP[b].JetThresholdEnergy / BHP[b].JetPower;
      jet_limit_myr = All.BHInternalTimestepFactor * bh_ffr_code_time_to_myr(tcode);
      if(!isfinite(jet_limit_myr) || !(jet_limit_myr > 0))
        terminate("BH_FFR: invalid jet timestep limit=%g Myr for ID=%llu", jet_limit_myr,
                  (unsigned long long)P[p].ID);
      if(jet_limit_myr < dt_limit_myr)
        dt_limit_myr = jet_limit_myr;
    }

  limited = bh_ffr_limit_integer_step_by_physical_myr(limited, dt_limit_myr);
  const integertime accuracy_limited = limited;

  /* If a channel still stores at least as many full thresholds as can be
   * released in one activation, request one bin finer than the timestep that
   * would otherwise be chosen by gas+internal-accuracy limits.  Anchor this
   * to accuracy_limited, not the particle's current bin: otherwise a persistent
   * backlog recursively halves the BH bin on every activation until the
   * minimum timeline step is reached. */
  int backlog = 0;
  if(BHP[b].WindThresholdEnergy > 0 &&
     BHP[b].WindEnergyBuffer >= All.BHMaxPacketsPerStep * BHP[b].WindThresholdEnergy)
    backlog = 1;
  if(BHP[b].JetThresholdEnergy > 0 &&
     BHP[b].JetEnergyBuffer >= All.BHMaxPacketsPerStep * BHP[b].JetThresholdEnergy)
    backlog = 1;

  integertime backlog_limit = accuracy_limited;
  int backlog_applied = 0;
  if(backlog && accuracy_limited > 2)
    {
      backlog_limit = accuracy_limited >> 1;
      if(backlog_limit < 2)
        backlog_limit = 2;

      if(backlog_limit < limited)
        {
          limited = backlog_limit;
          backlog_applied = 1;
        }
    }

  if(limited < raw_step)
    {
      const double wind_ratio =
          BHP[b].WindThresholdEnergy > 0 ? BHP[b].WindEnergyBuffer / BHP[b].WindThresholdEnergy : 0.0;
      const double jet_ratio =
          BHP[b].JetThresholdEnergy > 0 ? BHP[b].JetEnergyBuffer / BHP[b].JetThresholdEnergy : 0.0;

      printf("BH_FFR: timestep limit ID=%llu task=%d raw=%lld limited=%lld gasbin=%d fint=%g "
             "dtintMyr=%g dtwindMyr=%g dtjetMyr=%g accuracyLimited=%lld backlogLimit=%lld "
             "windBacklog=%g jetBacklog=%g backlog=%d backlogApplied=%d\n",
             (unsigned long long)P[p].ID, ThisTask, (long long)raw_step, (long long)limited, gas_bin,
             All.BHInternalTimestepFactor, reservoir_limit_myr, wind_limit_myr, jet_limit_myr,
             (long long)accuracy_limited, (long long)backlog_limit, wind_ratio, jet_ratio, backlog, backlog_applied);
      fflush(stdout);
    }

  return limited;
}

void bh_ffr_step(void)
{
  /* Particle array indices can change after domain/FoF reordering. Rebuild
   * this derived target cache from AREPO's current gravity-active list. */
  bh_ffr_build_active_list();

  /* Capture uses AREPO's generic MPI communication pattern. A BH may exist on
   * only one task, but every task in MPI_COMM_WORLD must enter the same
   * collectives in the same order. Therefore gate the collective phase on the
   * global, not local, number of active BHs. */
  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

  static int reservoir_self_test_done = 0;
  if(!reservoir_self_test_done)
    {
      bh_ffr_reservoir_self_test();
      bh_ffr_inner_self_test();
      bh_ffr_feedback_self_test();
      bh_ffr_jet_self_test();
      reservoir_self_test_done = 1;
    }

  bh_ffr_capture_resolved_gas();

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_step");

      if(P[p].Ti_Current != All.Ti_Current)
        terminate("BH_FFR: active particle ID=%llu is not drifted to Ti_Current=%lld", (unsigned long long)P[p].ID,
                  (long long)All.Ti_Current);

      const double dt_myr = bh_ffr_get_elapsed_time_myr(p);
      const double dt_code = bh_ffr_get_elapsed_time_code_time(p);
      if(BHP[b].LastProcessedTi != All.Ti_Current && (!(dt_myr > 0) || !(dt_code > 0)))
        terminate("BH_FFR: non-positive elapsed physical timestep for active particle ID=%llu", (unsigned long long)P[p].ID);

      /* Exact reservoir processing and Iteration-6 inner-flow partition are
       * committed together.  The update drains only the analytically
       * processable mass and assigns it exactly to horizon + wind channels. */
      bh_ffr_update_reservoir_state(p, dt_myr, dt_code);

      /* Persistent JetDir is a memory axis, not a physical spin vector.
       * Coherent reservoirs rotate it analytically toward DiscDir; incoherent
       * reservoirs leave it frozen. */
      bh_ffr_update_jet_direction(p, dt_myr);

      BHP[b].LastProcessedTi = All.Ti_Current;
    }

  /* Wind and jet channels use independent burst reservoirs and geometries.
   * They are applied sequentially so each channel solves its exact kinetic
   * packet against the gas state left by the preceding channel. */
  bh_ffr_inject_wind_feedback();
  bh_ffr_inject_jet_feedback();

  /* Iteration 11 merges synchronized BHs whose two proper accretion
   * apertures overlap. This happens after each object's local sub-grid
   * transaction/feedback for the current step. */
  bh_ffr_merge_close_black_holes();
}
