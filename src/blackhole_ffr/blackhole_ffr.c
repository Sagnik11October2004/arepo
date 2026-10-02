#include <math.h>
#include <stdint.h>
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

void bh_ffr_free_active_list(void)
{
  if(BHFFRActiveParticleList != NULL)
    myfree(BHFFRActiveParticleList);

  BHFFRActiveParticleList = NULL;
  NumActiveBHFFR = 0;
}

void bh_ffr_build_active_list(void)
{
  bh_ffr_free_active_list();

  if(TimeBinsGravity.NActiveParticles <= 0)
    return;

  BHFFRActiveParticleList =
      (int *)mymalloc("BHFFRActiveParticleList", TimeBinsGravity.NActiveParticles * sizeof(int));

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
    myfree(BHP);

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
      BHP = (struct bh_ffr_particle_data *)mymalloc("BHP", count * sizeof(struct bh_ffr_particle_data));
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

          if(!isfinite(BHP[b].BHMass) || BHP[b].BHMass < 0 || !isfinite(BHP[b].ReservoirMass) || BHP[b].ReservoirMass < 0 ||
             !isfinite(BHP[b].WindMassBuffer) || BHP[b].WindMassBuffer < 0 || !isfinite(BHP[b].WindMomentumBuffer) ||
             BHP[b].WindMomentumBuffer < 0 || !isfinite(BHP[b].WindEnergyBuffer) || BHP[b].WindEnergyBuffer < 0 ||
             !isfinite(BHP[b].JetEnergyBuffer) || BHP[b].JetEnergyBuffer < 0)
            terminate("BH_FFR: invalid mass/feedback buffer for particle ID=%llu in %s", (unsigned long long)P[i].ID, where);

          double jet_norm2 = 0;
          double disc_norm2 = 0;
          for(int k = 0; k < 3; k++)
            {
              if(!isfinite(BHP[b].Coherence[k]) || !isfinite(BHP[b].DiscDir[k]) || !isfinite(BHP[b].JetDir[k]))
                terminate("BH_FFR: non-finite direction state for particle ID=%llu in %s", (unsigned long long)P[i].ID, where);

              jet_norm2 += BHP[b].JetDir[k] * BHP[b].JetDir[k];
              disc_norm2 += BHP[b].DiscDir[k] * BHP[b].DiscDir[k];
            }

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

        BHP[b].ParticleID = P[i].ID;
        BHP[b].BHMass = P[i].Mass;
        BHP[b].ReservoirMass = 0;
        BHP[b].AccretionState = BH_FFR_STATE_UNINITIALIZED;

        bh_ffr_deterministic_axis(P[i].ID, BHP[b].DiscDir);
        for(int k = 0; k < 3; k++)
          BHP[b].JetDir[k] = BHP[b].DiscDir[k];

        b++;
      }

  bh_ffr_validate_state("initialization");

  long long local_count = NumBHFFR;
  long long global_count = 0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BH_FFR: initialized %lld Type-%d black-hole records; no accretion or feedback is active yet.\n", global_count,
             BH_FFR_PARTICLE_TYPE);
}
