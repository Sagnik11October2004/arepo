#include <math.h>

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

/*! Convert an AREPO integer-timeline interval to physical Myr.
 *
 * For cosmological runs the integer coordinate is dln(a), not physical time.
 * We reconstruct the endpoint scale factors exactly as predict.c does and use
 * AREPO's existing physical-time conversion helper. Non-cosmological runs use
 * the same helper with the corresponding code-time endpoints.
 */
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

/*! Physical elapsed time since this BH last executed the subgrid hook. */
double bh_ffr_get_elapsed_time_myr(int p)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_get_elapsed_time_myr");

  if(BHP[b].LastProcessedTi > All.Ti_Current)
    terminate("BH_FFR: LastProcessedTi=%lld exceeds Ti_Current=%lld for particle ID=%llu", (long long)BHP[b].LastProcessedTi,
              (long long)All.Ti_Current, (unsigned long long)P[p].ID);

  return bh_ffr_integer_interval_to_physical_myr(BHP[b].LastProcessedTi, All.Ti_Current);
}

/*! Store the minimum neighbouring gas hydro timebin returned by the gas pass.
 * A value of -1 means that no neighbour-derived limit is currently available.
 */
void bh_ffr_set_min_neighbour_timebin(int p, int timebin)
{
  const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_set_min_neighbour_timebin");

  if(timebin < -1 || timebin == 0 || timebin >= TIMEBINS)
    terminate("BH_FFR: invalid neighbour hydro timebin=%d for particle ID=%llu", timebin, (unsigned long long)P[p].ID);

  BHP[b].MinNeighbourHydroTimeBin = timebin;
}

/*! Apply the cached nearby-gas limit to a Type-5 gravity timestep.
 *
 * Iteration 2 deliberately does not invent a proxy for nearby gas activity.
 * MinNeighbourHydroTimeBin remains -1 until Iteration 3's distributed gas
 * search measures it. Once supplied, the normal AREPO gravity timebin
 * assignment automatically becomes equal to or finer than the fastest gas
 * neighbour inside the accretion aperture.
 */
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

/*! End-of-step FFR-MACER orchestration through Iteration 3.
 *
 * This advances only restart-safe synchronization bookkeeping. No gas mass,
 * BH mass, reservoir mass, momentum, radiation, wind, or jet state is changed.
 */
void bh_ffr_step(void)
{
  if(NumActiveBHFFR <= 0)
    return;

  /* Iteration 3: perform the read-only distributed gas-aperture pass before
   * advancing BH bookkeeping. This activates the nearby-gas timestep cache
   * without removing gas or changing any sub-grid mass/energy state. */
  bh_ffr_refresh_gas_neighbour_cache();

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = bh_ffr_compact_index_from_particle(p, "bh_ffr_step");

      if(P[p].Ti_Current != All.Ti_Current)
        terminate("BH_FFR: active particle ID=%llu is not drifted to Ti_Current=%lld", (unsigned long long)P[p].ID,
                  (long long)All.Ti_Current);

      const double dt_myr = bh_ffr_get_elapsed_time_myr(p);
      if(BHP[b].LastProcessedTi != All.Ti_Current && !(dt_myr > 0))
        terminate("BH_FFR: non-positive elapsed physical timestep for active particle ID=%llu", (unsigned long long)P[p].ID);

      BHP[b].LastProcessedTi = All.Ti_Current;
    }
}
