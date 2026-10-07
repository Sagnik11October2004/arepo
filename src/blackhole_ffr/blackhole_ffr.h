#ifndef BLACKHOLE_FFR_H
#define BLACKHOLE_FFR_H

#include "../main/allvars.h"

#define BH_FFR_PARTICLE_TYPE 5
/* Merger losers are temporarily retyped to a collisionless species so they
 * leave the live Type-5 set immediately, then compacted at the next normal
 * domain rearrangement. Snapshot I/O explicitly filters this exact
 * ID=0/Mass=0 tombstone signature. */
#define BH_FFR_MERGER_TOMBSTONE_TYPE 3
/* Version-1 cosmological DM selector. The default is Type 1, but
 * zoom ICs may place collisionless dark matter in several particle species.
 * Override this compile-time bitmask in Config.sh without changing All/BHP
 * restart layouts, e.g. BH_FFR_DM_TYPEMASK=2+4+8 for Types 1,2,3. */
#ifndef BH_FFR_DM_TYPEMASK
#define BH_FFR_DM_TYPEMASK (1 << 1)
#endif

/* Maximum nearest-neighbour payload carried by the SigmaDM gravity-tree
 * communication.  The model default needs O(40) particles, so 64 avoids the
 * previous 256-entry result structure on every exported BH while retaining
 * headroom for convergence tests.  Larger experiments can override this in
 * Config.sh without changing any persistent restart structure. */
#ifndef BH_FFR_DM_MAX_NEIGHBOURS
#define BH_FFR_DM_MAX_NEIGHBOURS 64
#endif

#if BH_FFR_DM_MAX_NEIGHBOURS < 2
#error "BH_FFR_DM_MAX_NEIGHBOURS must be at least 2"
#endif

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


/* Accretion-law selector for convergence benchmarks.  The default branch
 * behaviour is recovered with ACC_FFR + TARGET_RESERVOIR.  The shell-FFR
 * benchmark uses source-local constants in bh_ffr_capture.c so it does not
 * add Config.sh options or alter production parameter/restart layouts. */

enum bh_benchmark_accretion_model
{
  BH_BENCHMARK_ACC_TNG_BONDI = 0,
  BH_BENCHMARK_ACC_BOOSTED_BONDI = 1,
  BH_BENCHMARK_ACC_AM_BONDI = 2,
  BH_BENCHMARK_ACC_FFR = 3,
  BH_BENCHMARK_ACC_FFR_SHELL = 4,
  /* Convergence- and angular-momentum-corrected shell free-fall model.
   * The Weinberger-style shell m/t_ff estimator remains the backbone.
   * A coherent through-flow Mach factor is suppressed when the shell flow
   * is genuinely convergent, while a deliberately mild circularization
   * factor accounts for coherent shell angular momentum. */
  BH_BENCHMARK_ACC_CONVJ_SHELL_FFR = 5,
  BH_BENCHMARK_ACC_COUNT = 6
};

enum bh_benchmark_accretion_target
{
  BH_BENCHMARK_TARGET_RESERVOIR = 0,
  BH_BENCHMARK_TARGET_DIRECT = 1
};

enum bh_benchmark_boost_mode
{
  BH_BENCHMARK_BOOST_CONSTANT = 0,
  BH_BENCHMARK_BOOST_DENSITY = 1
};

/* Feedback is selected independently of the resolved accretion estimator.
 * NONE is useful for accretion-only convergence tests. TNG is the
 * Weinberger et al. thermal/kinetic law on the benchmark's controlled fixed
 * proper aperture. MACER selects the already-validated FFR-MACER wind+jet
 * implementation without changing its physics. */
enum bh_benchmark_feedback_model
{
  BH_BENCHMARK_FEEDBACK_NONE = 0,
  BH_BENCHMARK_FEEDBACK_TNG = 1,
  BH_BENCHMARK_FEEDBACK_MACER = 2
};

enum bh_benchmark_tng_feedback_mode
{
  BH_BENCHMARK_TNG_MODE_UNINITIALIZED = -1,
  BH_BENCHMARK_TNG_MODE_KINETIC = 0,
  BH_BENCHMARK_TNG_MODE_THERMAL = 1
};

struct bh_benchmark_environment
{
  double GasMass;
  double Density;
  double SoundSpeed;
  double RelativeSpeed;
  double Vphi;
  double HydrogenNumberDensity;
  double FFRRawRate;
  double FFRShellRawRate;
  double FFRShellGeometricRate;
  double FFRConvJShellRate;
  double ConvJFactor;
  double MachFactor;
  double AngularMomentumFactor;
  double ConvergenceFraction;
  double ShellBulkMach;
  double EffectiveMach;
  double CircularizationRatio;
  double CircularizationRadius;
  double ShellSoundSpeed;
  double ShellBulkSpeed;
  double ShellSpecificAngularMomentum;
};

struct bh_benchmark_rate_result
{
  double RawRate;
  double OperationalRate;
  double EddingtonRate;
  double BoostFactor;
  double AngularMomentumLimiter;
  double ConvJFactor;
  double MachFactor;
  double AngularMomentumFactor;
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

  /* Benchmark-only TNG feedback state.  The thermal buffer is normally
   * emptied on the same synchronization point; it remains persistent so the
   * active-cell-only public-AREPO coupling never discards energy.  The kinetic
   * buffer is the physical burst reservoir of the TNG model. */
  MyDouble BenchmarkMdotRaw;
  MyDouble BenchmarkMdotOperational;
  MyDouble BenchmarkMdotRealized;
  MyDouble BenchmarkSinkRateRatio;
  MyDouble BenchmarkCumulativeOperationalMass;
  MyDouble BenchmarkCumulativeRealizedMass;
  MyDouble BenchmarkCumulativeBHMassGrowth;

  MyDouble TNGThermalEnergyBuffer;
  MyDouble TNGKineticEnergyBuffer;
  MyDouble TNGFeedbackPower;
  MyDouble TNGKineticThresholdEnergy;
  MyDouble TNGEddingtonRatio;
  MyDouble TNGModeThreshold;
  MyDouble TNGCumulativeGeneratedEnergy;
  MyDouble TNGCumulativeInjectedEnergy;
  MyDouble TNGActiveTargetMassFraction;
  MyDouble TNGThermalBufferAgeCodeTime;
  MyDouble TNGKineticBufferAgeCodeTime;

  int AccretionState;
  int TNGFeedbackMode;
  int TNGThermalBufferAgeTransactions;
  int TNGKineticBufferAgeTransactions;
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

/* Small fixed-capacity container for deterministic discrete adaptive-radius
 * candidates.  Iteration A only provides validated infrastructure; no live
 * accretion or feedback path selects from this grid yet.  Radii are proper
 * code lengths, matching BHAccretionRadius/BHFeedbackRadius semantics. */
#define BH_FFR_MAX_ADAPTIVE_RADII 8
struct bh_ffr_discrete_radius_grid
{
  int Count;
  MyDouble Radius[BH_FFR_MAX_ADAPTIVE_RADII];
};

/*! Restart record for the transient dynamical-friction environment.
 *
 * This remains separate from bh_ffr_particle_data so the DM cache is not
 * replicated inside every persistent BH record.  The benchmark branch now
 * versions its extended BHP layout explicitly as native restart version 4.
 * RestartFlag=1 stores only cache entries owned by the local rank, then
 * reconstructs the replicated ID-addressable cache after all rank-local
 * restart files have been read.
 */
struct bh_ffr_df_restart_entry
{
  MyIDType ID;
  integertime Ti;
  MyDouble Pos[3];
  double RhoDM;
  double SigmaDM;
  double MeanVel[3];
  int Valid;
};

void bh_ffr_allocate_state(int count);
void bh_ffr_free_state(void);
void bh_ffr_initialize_particles(void);
void bh_ffr_rebuild_state_after_particle_changes(void);
void bh_ffr_validate_state(const char *where);

void bh_ffr_build_discrete_radius_grid(double base_radius, const double *factors, int count,
                                       struct bh_ffr_discrete_radius_grid *grid);
int bh_ffr_select_smallest_resolved_radius(const struct bh_ffr_discrete_radius_grid *grid,
                                           const long long *counts, long long min_count,
                                           int *underresolved);
void bh_ffr_adaptive_radius_self_test(void);

void bh_ffr_build_active_list(void);
void bh_ffr_free_active_list(void);

void bh_ffr_step(void);
void bh_ffr_collect_gas_environment(struct bh_ffr_gas_search_result *results);
void bh_ffr_refresh_gas_neighbour_cache(void);
void bh_ffr_prepare_dm_environment_search(void);
void bh_ffr_apply_cached_dynamical_friction(int p, double dt_code);
double bh_ffr_get_cached_dynamical_friction_timescale_code(int p);
int bh_ffr_df_restart_local_count(void);
void bh_ffr_df_restart_export_local(struct bh_ffr_df_restart_entry *out, int count);
void bh_ffr_df_restart_import_local(const struct bh_ffr_df_restart_entry *in, int count);
void bh_ffr_df_restart_replicate(void);
void bh_ffr_merge_close_black_holes(void);
void bh_ffr_capture_resolved_gas(void);
void bh_ffr_capture_self_test(void);
const char *bh_benchmark_accretion_model_name(int model);
double bh_benchmark_eddington_rate_code(double bh_mass);
void bh_benchmark_compute_accretion(const struct bh_benchmark_environment *env, double bh_mass,
                                    struct bh_benchmark_rate_result *out);
void bh_benchmark_compute_all_raw_rates(const struct bh_benchmark_environment *env, double bh_mass,
                                        double raw_rates[BH_BENCHMARK_ACC_COUNT],
                                        double *boost_factor, double *am_limiter,
                                        double *env_factor);
void bh_benchmark_accretion_self_test(void);
void bh_benchmark_tng_feedback_accumulate(void);
void bh_benchmark_tng_feedback_inject(void);
void bh_benchmark_tng_feedback_self_test(void);
double bh_ffr_central_mass_code(int p);
double bh_ffr_reservoir_timescale_myr(double bh_mass, double reservoir_mass);
double bh_ffr_reservoir_processable_mass(double reservoir_mass, double dt_myr, double tau_myr);
double bh_ffr_eddington_rate_code(double bh_mass);
double bh_ffr_truncation_to_hot_radius_ratio(double processed_edd_ratio);
int bh_ffr_classify_accretion_state(int previous_state, double processed_edd_ratio, double rtr_over_rhot);
void bh_ffr_update_reservoir_state(int p, double dt_myr, double dt_code, double disk_time_myr);
void bh_ffr_reservoir_self_test(void);
void bh_ffr_apply_inner_flow(int p, double processable_mass, double dt_code);
void bh_ffr_inner_self_test(void);
void bh_ffr_inject_wind_feedback(void);
void bh_ffr_feedback_self_test(void);
void bh_ffr_update_jet_direction(int p, double dt_myr, double disk_time_myr);
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