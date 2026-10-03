#ifndef BLACKHOLE_FFR_H
#define BLACKHOLE_FFR_H

#include "../main/allvars.h"

#define BH_FFR_PARTICLE_TYPE 5
/* Version-1 cosmological DM selector. The current AREPO test/production
 * configuration uses Type-1 as dark matter. Keep this compile-time rather
 * than adding an All/BHP field so native restart layouts remain unchanged. */
#define BH_FFR_DM_PARTICLE_TYPE 1

/* Iteration-5 state-machine closures from the design specification. The
 * nominal physical boundaries are dotm=0.02 and Rtr/Rhot=1; the wider
 * entry/exit values provide hysteresis against state chatter. */
#define BH_FFR_COLD_NOMINAL_EDD_RATIO 0.020
#define BH_FFR_COLD_ENTER_EDD_RATIO 0.022
#define BH_FFR_COLD_EXIT_EDD_RATIO 0.018
#define BH_FFR_ADIOS_ENTER_RTR_FACTOR 1.20
#define BH_FFR_ADIOS_EXIT_RTR_FACTOR 0.80

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
  /* Exact finite-step reservoir processing rate. Iteration 6 partitions
   * this rate conservatively into MdotHorizon + MdotWind. */
  MyDouble MdotProcessed;
  MyDouble MdotEddington;
  MyDouble ProcessedEddingtonRatio;
  /* Reserved Iteration-5 ABI slot. Cold blending is deferred to later
   * feedback iterations; retain this member so native restart files written
   * by the first Iteration-5 implementation keep the same BHP layout. */
  MyDouble ColdBlendWeight;
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
void bh_ffr_rebuild_state_after_particle_changes(void);
void bh_ffr_validate_state(const char *where);

void bh_ffr_build_active_list(void);
void bh_ffr_free_active_list(void);

void bh_ffr_step(void);
void bh_ffr_collect_gas_environment(struct bh_ffr_gas_search_result *results);
void bh_ffr_refresh_gas_neighbour_cache(void);
void bh_ffr_prepare_dm_environment_search(void);
void bh_ffr_apply_cached_dynamical_friction(int p, double dt_code);
void bh_ffr_merge_close_black_holes(void);
void bh_ffr_capture_resolved_gas(void);
void bh_ffr_capture_self_test(void);
double bh_ffr_reservoir_processable_mass(double reservoir_mass, double dt_myr, double tau_myr);
double bh_ffr_eddington_rate_code(double bh_mass);
double bh_ffr_truncation_to_hot_radius_ratio(double processed_edd_ratio);
int bh_ffr_classify_accretion_state(int previous_state, double processed_edd_ratio, double rtr_over_rhot);
void bh_ffr_update_reservoir_state(int p, double dt_myr, double dt_code);
void bh_ffr_reservoir_self_test(void);
void bh_ffr_apply_inner_flow(int p, double processable_mass, double dt_code);
void bh_ffr_inner_self_test(void);
void bh_ffr_inject_wind_feedback(void);
void bh_ffr_feedback_self_test(void);
void bh_ffr_update_jet_direction(int p, double dt_myr);
void bh_ffr_inject_jet_feedback(void);
void bh_ffr_jet_self_test(void);
int bh_ffr_seed_from_fof(void);
double bh_ffr_integer_interval_to_physical_myr(integertime ti0, integertime ti1);
double bh_ffr_integer_interval_to_physical_code_time(integertime ti0, integertime ti1);
double bh_ffr_get_elapsed_time_myr(int p);
double bh_ffr_get_elapsed_time_code_time(int p);
integertime bh_ffr_limit_gravity_timestep(int p, integertime ti_step);
void bh_ffr_set_min_neighbour_timebin(int p, int timebin);

void bh_ffr_domain_exchange_begin(struct bh_ffr_domain_exchange_context *ctx);
void bh_ffr_domain_exchange_finish(struct bh_ffr_domain_exchange_context *ctx);

#endif /* BLACKHOLE_FFR_H */
