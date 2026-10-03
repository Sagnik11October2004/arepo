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

integertime bh_ffr_limit_gravity_timestep(int p, integertime ti_step)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_limit_gravity_timestep");
  const int gas_bin = BHP[b].MinNeighbourHydroTimeBin;

  if(gas_bin < 0)
    return ti_step;

  if(gas_bin == 0 || gas_bin >= TIMEBINS)
    terminate("BH_FFR: corrupt neighbour hydro timebin=%d for particle ID=%llu", gas_bin, (unsigned long long)P[p].ID);

  const integertime gas_ti_step = ((integertime)1) << gas_bin;
  return gas_ti_step < ti_step ? gas_ti_step : ti_step;
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

      /* Iteration 5 computes the exact finite-step processable reservoir mass
       * and state diagnostics, but deliberately does not drain ReservoirMass.
       * Iteration 6 will atomically partition this candidate into horizon and
       * wind mass before committing the reservoir transaction. */
      bh_ffr_update_reservoir_state(p, dt_myr, dt_code);

      BHP[b].LastProcessedTi = All.Ti_Current;
    }
}
