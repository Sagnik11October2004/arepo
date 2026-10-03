#ifndef BLACKHOLE_FFR_H
#define BLACKHOLE_FFR_H

#include "../main/allvars.h"

#define BH_FFR_PARTICLE_TYPE 5

enum bh_ffr_accretion_state
{
  BH_FFR_STATE_UNINITIALIZED = -1,
  BH_FFR_STATE_ADIOS = 0,
  BH_FFR_STATE_TRUNCATED = 1,
  BH_FFR_STATE_COLD = 2
};

/*! Persistent sub-grid state for one Type-5 FFR-MACER black hole.
 *
 * The coherence vector stores mass-weighted orientation memory, not a
 * physical angular-momentum magnitude. WindMomentumBuffer is a scalar
 * unresolved impulse budget; a later bipolar event will split it into
 * equal and opposite lobes.
 */
struct bh_ffr_particle_data
{
  MyIDType ParticleID;

  /* Restart-safe synchronization bookkeeping. LastProcessedTi prevents a
   * spurious full timebin update on startup/restart. The neighbour timebin
   * is populated by the later gas tree pass and limits the next BH step. */
  integertime LastProcessedTi;
  int MinNeighbourHydroTimeBin;

  MyDouble BHMass;
  MyDouble ReservoirMass;

  MyDouble Coherence[3];
  MyDouble DiscDir[3];
  MyDouble JetDir[3];

  MyDouble WindMassBuffer;
  MyDouble WindMomentumBuffer;
  MyDouble WindEnergyBuffer;
  MyDouble JetEnergyBuffer;

  MyDouble MdotSupply;
  MyDouble MdotProcessed;
  MyDouble MdotHorizon;
  MyDouble MdotWind;
  MyDouble BolometricLuminosity;
  MyDouble WindPower;
  MyDouble JetPower;
  MyDouble SigmaDM;
  MyDouble WindThresholdEnergy;
  MyDouble JetThresholdEnergy;

  int AccretionState;
};

extern struct bh_ffr_particle_data *BHP;
extern int NumBHFFR;

extern int *BHFFRActiveParticleList;
extern int NumActiveBHFFR;

/*! Read-only distributed gas-aperture result for one active BH.
 *
 * Counts and masses are accumulated independently inside the fixed proper
 * accretion and feedback apertures. MinHydroTimeBin refers only to gas inside
 * the accretion aperture; -1 means no finite neighbour-derived limit.
 */
struct bh_ffr_gas_search_result
{
  MyDouble AccretionMass;
  MyDouble FeedbackMass;
  long long AccretionCount;
  long long FeedbackCount;
  int MinHydroTimeBin;
};

struct bh_ffr_domain_exchange_context
{
  struct bh_ffr_particle_data *Received;
  int NumReceived;
};

void bh_ffr_allocate_state(int count);
void bh_ffr_free_state(void);
void bh_ffr_initialize_particles(void);
void bh_ffr_validate_state(const char *where);

void bh_ffr_build_active_list(void);
void bh_ffr_free_active_list(void);

void bh_ffr_step(void);
void bh_ffr_collect_gas_environment(struct bh_ffr_gas_search_result *results);
void bh_ffr_refresh_gas_neighbour_cache(void);
void bh_ffr_prepare_dm_environment_search(void);
double bh_ffr_integer_interval_to_physical_myr(integertime ti0, integertime ti1);
double bh_ffr_get_elapsed_time_myr(int p);
integertime bh_ffr_limit_gravity_timestep(int p, integertime ti_step);
void bh_ffr_set_min_neighbour_timebin(int p, int timebin);

void bh_ffr_domain_exchange_begin(struct bh_ffr_domain_exchange_context *ctx);
void bh_ffr_domain_exchange_finish(struct bh_ffr_domain_exchange_context *ctx);

#endif /* BLACKHOLE_FFR_H */
