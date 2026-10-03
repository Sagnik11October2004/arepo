#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

struct bh_ffr_particle_data *BHP = NULL;
int NumBHFFR = 0;

int *BHFFRActiveParticleList = NULL;
int NumActiveBHFFR = 0;

static uint64_t bh_ffr_hash64(uint64_t x)
{
  x += UINT64_C(0x9e3779b97f4a7c15);
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  return x ^ (x >> 31);
}

static double bh_ffr_hash_unit_interval(uint64_t x)
{
  return (double)(bh_ffr_hash64(x) >> 11) / 9007199254740992.0;
}

static void bh_ffr_deterministic_axis(MyIDType id, MyDouble axis[3])
{
  const uint64_t key = (uint64_t)id;
  const double z = 2.0 * bh_ffr_hash_unit_interval(key) - 1.0;
  const double phi = 6.2831853071795864769 * bh_ffr_hash_unit_interval(key ^ UINT64_C(0xd1b54a32d192ed03));
  const double r = sqrt(dmax(0.0, 1.0 - z * z));

  axis[0] = r * cos(phi);
  axis[1] = r * sin(phi);
  axis[2] = z;
}

static void bh_ffr_initialize_record(struct bh_ffr_particle_data *bh, int p, integertime last_processed_ti)
{
  memset(bh, 0, sizeof(*bh));

  bh->ParticleID = P[p].ID;
  bh->LastProcessedTi = last_processed_ti;
  bh->MinNeighbourHydroTimeBin = -1;
  bh->BHMass = P[p].Mass;
  bh->ReservoirMass = 0;
  bh->AccretionState = BH_FFR_STATE_UNINITIALIZED;

  bh_ffr_deterministic_axis(P[p].ID, bh->DiscDir);
  for(int k = 0; k < 3; k++)
    bh->JetDir[k] = bh->DiscDir[k];
}

void bh_ffr_free_active_list(void)
{
  if(BHFFRActiveParticleList != NULL)
    free(BHFFRActiveParticleList);

  BHFFRActiveParticleList = NULL;
  NumActiveBHFFR = 0;
}

void bh_ffr_build_active_list(void)
{
  bh_ffr_free_active_list();

  if(TimeBinsGravity.NActiveParticles <= 0)
    return;

  /* The active list can remain live while FoF seeding reconstructs timebins
   * and tears down/rebuilds the gas neighbour tree. It therefore must not
   * occupy AREPO's LIFO mymalloc arena behind movable tree blocks. */
  BHFFRActiveParticleList =
      (int *)malloc(TimeBinsGravity.NActiveParticles * sizeof(int));
  if(BHFFRActiveParticleList == NULL)
    terminate("BH_FFR: failed to allocate active-particle list for %d gravity-active particles",
              TimeBinsGravity.NActiveParticles);

  for(int n = 0; n < TimeBinsGravity.NActiveParticles; n++)
    {
      const int p = TimeBinsGravity.ActiveParticleList[n];

      if(p < 0)
        continue;
      if(p >= NumPart)
        terminate("BH_FFR: active gravity-list index %d is outside NumPart=%d", p, NumPart);

      if(P[p].Type == BH_FFR_PARTICLE_TYPE)
        {
          const int b = P[p].BHDataIndex;
          if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
            terminate("BH_FFR: invalid compact state for active Type-5 particle ID=%llu", (unsigned long long)P[p].ID);

          BHFFRActiveParticleList[NumActiveBHFFR++] = p;
        }
    }

  if(NumActiveBHFFR == 0)
    bh_ffr_free_active_list();
}

void bh_ffr_free_state(void)
{
  bh_ffr_free_active_list();

  if(BHP != NULL)
    free(BHP);

  BHP = NULL;
  NumBHFFR = 0;
}

void bh_ffr_allocate_state(int count)
{
  if(count < 0)
    terminate("BH_FFR: negative compact-state size %d", count);

  bh_ffr_free_state();

  if(count > 0)
    {
      /* BHP must survive domain decomposition, whose AREPO allocator stack
       * contains unrelated non-movable blocks. Keep this persistent auxiliary
       * state on the C heap rather than inside the mymalloc LIFO arena. */
      BHP = (struct bh_ffr_particle_data *)malloc(count * sizeof(struct bh_ffr_particle_data));
      if(BHP == NULL)
        terminate("BH_FFR: failed to allocate compact state for %d black holes", count);
      memset(BHP, 0, count * sizeof(struct bh_ffr_particle_data));
    }

  NumBHFFR = count;
}

void bh_ffr_validate_state(const char *where)
{
  int local_bh_count = 0;
  unsigned char *seen = NULL;

  if(NumBHFFR > 0)
    {
      seen = (unsigned char *)mymalloc("BHFFRSeen", NumBHFFR * sizeof(unsigned char));
      memset(seen, 0, NumBHFFR * sizeof(unsigned char));
    }

  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type == BH_FFR_PARTICLE_TYPE)
        {
          local_bh_count++;

          const int b = P[i].BHDataIndex;
          if(b < 0 || b >= NumBHFFR)
            terminate("BH_FFR: invalid BHDataIndex=%d for particle ID=%llu in %s", b, (unsigned long long)P[i].ID, where);

          if(seen[b])
            terminate("BH_FFR: duplicate compact index %d in %s", b, where);
          seen[b] = 1;

          if(BHP[b].ParticleID != P[i].ID)
            terminate("BH_FFR: ID/index mismatch for particle ID=%llu in %s", (unsigned long long)P[i].ID, where);

          if(BHP[b].LastProcessedTi < 0 || BHP[b].LastProcessedTi > All.Ti_Current)
            terminate("BH_FFR: invalid LastProcessedTi=%lld at Ti_Current=%lld for particle ID=%llu in %s",
                      (long long)BHP[b].LastProcessedTi, (long long)All.Ti_Current, (unsigned long long)P[i].ID, where);

          if(BHP[b].MinNeighbourHydroTimeBin < -1 || BHP[b].MinNeighbourHydroTimeBin == 0 ||
             BHP[b].MinNeighbourHydroTimeBin >= TIMEBINS)
            terminate("BH_FFR: invalid minimum neighbour hydro timebin=%d for particle ID=%llu in %s",
                      BHP[b].MinNeighbourHydroTimeBin, (unsigned long long)P[i].ID, where);

          if(!isfinite(BHP[b].BHMass) || BHP[b].BHMass < 0 || !isfinite(BHP[b].ReservoirMass) || BHP[b].ReservoirMass < 0 ||
             !isfinite(BHP[b].WindMassBuffer) || BHP[b].WindMassBuffer < 0 || !isfinite(BHP[b].WindMomentumBuffer) ||
             BHP[b].WindMomentumBuffer < 0 || !isfinite(BHP[b].WindEnergyBuffer) || BHP[b].WindEnergyBuffer < 0 ||
             !isfinite(BHP[b].JetEnergyBuffer) || BHP[b].JetEnergyBuffer < 0 || !isfinite(BHP[b].MdotSupply) ||
             BHP[b].MdotSupply < 0 || !isfinite(BHP[b].MdotProcessed) || BHP[b].MdotProcessed < 0 ||
             !isfinite(BHP[b].MdotEddington) || BHP[b].MdotEddington < 0 ||
             !isfinite(BHP[b].ProcessedEddingtonRatio) || BHP[b].ProcessedEddingtonRatio < 0 ||
             !isfinite(BHP[b].MdotHorizon) || BHP[b].MdotHorizon < 0 ||
             !isfinite(BHP[b].MdotWind) || BHP[b].MdotWind < 0 ||
             !isfinite(BHP[b].BolometricLuminosity) || BHP[b].BolometricLuminosity < 0 ||
             !isfinite(BHP[b].WindPower) || BHP[b].WindPower < 0 ||
             !isfinite(BHP[b].JetPower) || BHP[b].JetPower < 0 ||
             !isfinite(BHP[b].WindThresholdEnergy) || BHP[b].WindThresholdEnergy < 0 ||
             !isfinite(BHP[b].JetThresholdEnergy) || BHP[b].JetThresholdEnergy < 0 ||
             !isfinite(BHP[b].ColdBlendWeight) || BHP[b].ColdBlendWeight < 0 || BHP[b].ColdBlendWeight > 1)
            terminate("BH_FFR: invalid mass/rate/energetics buffer for particle ID=%llu in %s",
                      (unsigned long long)P[i].ID, where);

          if(BHP[b].AccretionState < BH_FFR_STATE_UNINITIALIZED || BHP[b].AccretionState > BH_FFR_STATE_COLD)
            terminate("BH_FFR: invalid accretion state=%d for particle ID=%llu in %s", BHP[b].AccretionState,
                      (unsigned long long)P[i].ID, where);

          /* Old Iteration-5 native restarts can legitimately contain a
           * diagnostic MdotProcessed with zero horizon/wind rates.  Enforce
           * the Iteration-6 partition once either inner channel is populated;
           * the live transaction itself checks exact closure unconditionally. */
          if(BHP[b].MdotHorizon > 0 || BHP[b].MdotWind > 0)
            {
              const double inner_rate_scale = fmax(1.0, BHP[b].MdotProcessed);
              if(fabs(BHP[b].MdotProcessed - (BHP[b].MdotHorizon + BHP[b].MdotWind)) > 2.0e-10 * inner_rate_scale)
                terminate("BH_FFR: inner-flow rate partition mismatch for particle ID=%llu in %s: proc=%g H=%g wind=%g",
                          (unsigned long long)P[i].ID, where, BHP[b].MdotProcessed, BHP[b].MdotHorizon, BHP[b].MdotWind);
            }

          const double expected_dyn_mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
          const double dyn_scale = dmax(1.0, dmax(fabs(expected_dyn_mass), fabs(P[i].Mass)));
          if(!isfinite(P[i].Mass) || P[i].Mass < 0 || fabs(P[i].Mass - expected_dyn_mass) > 1.0e-10 * dyn_scale)
            terminate("BH_FFR: dynamical-mass mismatch for particle ID=%llu in %s: P.Mass=%g components=%g",
                      (unsigned long long)P[i].ID, where, P[i].Mass, expected_dyn_mass);

          double jet_norm2 = 0;
          double disc_norm2 = 0;
          double coherence_norm2 = 0;
          for(int k = 0; k < 3; k++)
            {
              if(!isfinite(BHP[b].Coherence[k]) || !isfinite(BHP[b].DiscDir[k]) || !isfinite(BHP[b].JetDir[k]))
                terminate("BH_FFR: non-finite direction state for particle ID=%llu in %s", (unsigned long long)P[i].ID, where);

              jet_norm2 += BHP[b].JetDir[k] * BHP[b].JetDir[k];
              disc_norm2 += BHP[b].DiscDir[k] * BHP[b].DiscDir[k];
              coherence_norm2 += BHP[b].Coherence[k] * BHP[b].Coherence[k];
            }

          const double coherence_norm = sqrt(coherence_norm2);
          if(coherence_norm > BHP[b].ReservoirMass * (1.0 + 1.0e-10))
            terminate("BH_FFR: coherence norm=%g exceeds reservoir mass=%g for particle ID=%llu in %s", coherence_norm,
                      BHP[b].ReservoirMass, (unsigned long long)P[i].ID, where);

          if(fabs(jet_norm2 - 1.0) > 1.0e-5 || fabs(disc_norm2 - 1.0) > 1.0e-5)
            terminate("BH_FFR: non-unit persistent axis for particle ID=%llu in %s", (unsigned long long)P[i].ID, where);
        }
      else if(P[i].BHDataIndex != -1)
        terminate("BH_FFR: non-BH particle ID=%llu carries BHDataIndex=%d in %s", (unsigned long long)P[i].ID, P[i].BHDataIndex,
                  where);
    }

  if(local_bh_count != NumBHFFR)
    terminate("BH_FFR: local Type-5 count=%d differs from NumBHFFR=%d in %s", local_bh_count, NumBHFFR, where);

  if(seen != NULL)
    myfree(seen);
}

void bh_ffr_initialize_particles(void)
{
  int count = 0;

  for(int i = 0; i < NumPart; i++)
    {
      P[i].BHDataIndex = -1;
      if(P[i].Type == BH_FFR_PARTICLE_TYPE)
        count++;
    }

  bh_ffr_allocate_state(count);

  int b = 0;
  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE)
      {
        P[i].BHDataIndex = b;

        bh_ffr_initialize_record(&BHP[b], i, 0);
        b++;
      }

  bh_ffr_validate_state("initialization");

  long long local_count = NumBHFFR;
  long long global_count = 0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BH_FFR: initialized %lld Type-%d black-hole records.\n", global_count, BH_FFR_PARTICLE_TYPE);
}


void bh_ffr_rebuild_state_after_particle_changes(void)
{
  bh_ffr_free_active_list();

  struct bh_ffr_particle_data *old_data = BHP;
  const int old_count = NumBHFFR;

  int new_count = 0;
  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE)
      new_count++;

  struct bh_ffr_particle_data *new_data = NULL;
  unsigned char *old_used = NULL;

  if(new_count > 0)
    {
      new_data = (struct bh_ffr_particle_data *)malloc(new_count * sizeof(*new_data));
      if(new_data == NULL)
        terminate("BH_FFR: failed to allocate rebuilt compact state for %d black holes", new_count);
    }

  if(old_count > 0)
    {
      old_used = (unsigned char *)mymalloc("BHFFRRebuildOldUsed", old_count * sizeof(*old_used));
      memset(old_used, 0, old_count * sizeof(*old_used));
    }

  int bnew = 0;
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type != BH_FFR_PARTICLE_TYPE)
        {
          P[i].BHDataIndex = -1;
          continue;
        }

      int found = -1;
      const int old_index = P[i].BHDataIndex;

      if(old_index >= 0 && old_index < old_count && !old_used[old_index] && old_data[old_index].ParticleID == P[i].ID)
        found = old_index;
      else
        for(int b = 0; b < old_count; b++)
          if(!old_used[b] && old_data[b].ParticleID == P[i].ID)
            {
              found = b;
              break;
            }

      if(found >= 0)
        {
          new_data[bnew] = old_data[found];
          old_used[found] = 1;
        }
      else
        bh_ffr_initialize_record(&new_data[bnew], i, All.Ti_Current);

      P[i].BHDataIndex = bnew;
      bnew++;
    }

  if(old_used != NULL)
    myfree(old_used);
  if(old_data != NULL)
    free(old_data);

  BHP = new_data;
  NumBHFFR = new_count;

  bh_ffr_validate_state("particle-change rebuild");
}
