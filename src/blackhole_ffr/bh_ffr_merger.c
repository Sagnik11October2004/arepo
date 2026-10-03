#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Iteration 11: simple unresolved BH mergers.
 *
 * Only gravity-active BHs at the current synchronization point participate.
 * Two BHs are linked when their proper accretion apertures overlap:
 *
 *      d < 2 R_acc .
 *
 * Every connected overlap component is collapsed deterministically onto the
 * lowest particle ID. The survivor keeps its position, while its velocity is
 * the exact dynamical-mass-weighted mean. Persistent unresolved mass/energy
 * reservoirs are summed. The loser is marked with AREPO's standard deleted
 * particle convention (ID=0, Mass=0) and removed from the gravity timebin.
 *
 * This is intentionally a simple sub-grid coalescence rule: no binding-energy,
 * relative-velocity, hardening, recoil, or gravitational-wave delay criterion
 * is imposed in version 1.
 */

struct bh_ffr_merge_summary
{
  MyIDType ID;
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyDouble DynMass;
  struct bh_ffr_particle_data State;
};

static int bh_ffr_merge_compare_id(const void *aa, const void *bb)
{
  const struct bh_ffr_merge_summary *a = aa;
  const struct bh_ffr_merge_summary *b = bb;
  if(a->ID < b->ID)
    return -1;
  if(a->ID > b->ID)
    return 1;
  return 0;
}

static int bh_ffr_merge_find(int *parent, int x)
{
  int r = x;
  while(parent[r] != r)
    r = parent[r];

  while(parent[x] != x)
    {
      const int n = parent[x];
      parent[x] = r;
      x = n;
    }

  return r;
}

static int bh_ffr_merge_find_const(const int *parent, int x)
{
  while(parent[x] != x)
    x = parent[x];
  return x;
}

static void bh_ffr_merge_union(int *parent, int a, int b)
{
  a = bh_ffr_merge_find(parent, a);
  b = bh_ffr_merge_find(parent, b);
  if(a == b)
    return;

  /* Summaries are sorted by particle ID; attaching the larger root to the
   * smaller root makes the lowest ID the deterministic survivor. */
  if(a < b)
    parent[b] = a;
  else
    parent[a] = b;
}

static void bh_ffr_normalize_axis(double v[3], const double fallback[3])
{
  double n2 = 0.0;
  for(int k = 0; k < 3; k++)
    n2 += v[k] * v[k];

  if(n2 <= 1.0e-30)
    {
      for(int k = 0; k < 3; k++)
        v[k] = fallback[k];
      return;
    }

  const double inv = 1.0 / sqrt(n2);
  for(int k = 0; k < 3; k++)
    v[k] *= inv;
}

static int bh_ffr_local_particle_from_id(MyIDType id)
{
  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE && P[i].ID == id)
      return i;
  return -1;
}

static void bh_ffr_build_merged_state(const struct bh_ffr_merge_summary *all, const int *parent, int root, int nall,
                                      struct bh_ffr_particle_data *merged, double vel[3], int *nmembers)
{
  memset(merged, 0, sizeof(*merged));
  *nmembers = 0;

  const struct bh_ffr_merge_summary *survivor = &all[root];
  merged->ParticleID = survivor->ID;
  merged->LastProcessedTi = All.Ti_Current;
  merged->MinNeighbourHydroTimeBin = -1;
  merged->AccretionState = BH_FFR_STATE_UNINITIALIZED;

  double dynmass = 0.0;
  double jetvec[3] = {0, 0, 0};
  double sigma_weight = 0.0;

  for(int i = 0; i < nall; i++)
    if(bh_ffr_merge_find_const(parent, i) == root)
      {
        const struct bh_ffr_particle_data *s = &all[i].State;
        const double w = all[i].DynMass;

        (*nmembers)++;
        dynmass += w;

        merged->BHMass += s->BHMass;
        merged->ReservoirMass += s->ReservoirMass;
        merged->WindMassBuffer += s->WindMassBuffer;
        merged->WindMomentumBuffer += s->WindMomentumBuffer;
        merged->WindEnergyBuffer += s->WindEnergyBuffer;
        merged->JetEnergyBuffer += s->JetEnergyBuffer;

        for(int k = 0; k < 3; k++)
          {
            merged->Coherence[k] += s->Coherence[k];
            vel[k] += w * all[i].Vel[k];
            jetvec[k] += s->BHMass * s->JetDir[k];
          }

        sigma_weight += w * s->SigmaDM;
      }

  if(!(dynmass > 0) || !isfinite(dynmass))
    terminate("BH_FFR: invalid merged dynamical mass=%g", dynmass);

  for(int k = 0; k < 3; k++)
    vel[k] /= dynmass;

  const double ledger = merged->BHMass + merged->ReservoirMass + merged->WindMassBuffer;
  if(fabs(ledger - dynmass) > 2.0e-10 * fmax(1.0, dynmass))
    terminate("BH_FFR: merger mass ledger mismatch components=%g particle-sum=%g", ledger, dynmass);

  double survivor_disc[3], survivor_jet[3];
  for(int k = 0; k < 3; k++)
    {
      survivor_disc[k] = survivor->State.DiscDir[k];
      survivor_jet[k] = survivor->State.JetDir[k];
      merged->DiscDir[k] = merged->Coherence[k];
      merged->JetDir[k] = jetvec[k];
    }

  bh_ffr_normalize_axis(merged->DiscDir, survivor_disc);
  bh_ffr_normalize_axis(merged->JetDir, survivor_jet);

  merged->SigmaDM = sigma_weight / dynmass;
  merged->MdotEddington = bh_ffr_eddington_rate_code(merged->BHMass);

  /* Instantaneous rates/powers are reset because they refer to the two
   * pre-merger objects. Conserved unresolved mass/energy buffers above are
   * retained exactly and will continue from the merged object next step. */
  merged->MdotSupply = 0.0;
  merged->MdotProcessed = 0.0;
  merged->ProcessedEddingtonRatio = 0.0;
  merged->ColdBlendWeight = 0.0;
  merged->MdotHorizon = 0.0;
  merged->MdotWind = 0.0;
  merged->BolometricLuminosity = 0.0;
  merged->WindPower = 0.0;
  merged->JetPower = 0.0;
  merged->WindThresholdEnergy = 0.0;
  merged->JetThresholdEnergy = 0.0;
}

void bh_ffr_merge_close_black_holes(void)
{
  const int local_n = NumActiveBHFFR;

  int *counts = malloc(NTask * sizeof(int));
  int *displs = malloc(NTask * sizeof(int));
  int *byte_counts = malloc(NTask * sizeof(int));
  int *byte_displs = malloc(NTask * sizeof(int));
  if(counts == NULL || displs == NULL || byte_counts == NULL || byte_displs == NULL)
    terminate("BH_FFR: failed to allocate merger MPI counts");

  MPI_Allgather(&local_n, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  int nall = 0;
  for(int task = 0; task < NTask; task++)
    {
      displs[task] = nall;
      nall += counts[task];
    }

  if(nall < 2)
    {
      free(byte_displs);
      free(byte_counts);
      free(displs);
      free(counts);
      return;
    }

  struct bh_ffr_merge_summary *local =
      malloc((local_n > 0 ? local_n : 1) * sizeof(struct bh_ffr_merge_summary));
  struct bh_ffr_merge_summary *all =
      malloc(nall * sizeof(struct bh_ffr_merge_summary));
  if(local == NULL || all == NULL)
    terminate("BH_FFR: failed to allocate merger summaries");

  for(int n = 0; n < local_n; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = P[p].BHDataIndex;

      local[n].ID = P[p].ID;
      local[n].DynMass = P[p].Mass;
      for(int k = 0; k < 3; k++)
        {
          local[n].Pos[k] = P[p].Pos[k];
          local[n].Vel[k] = P[p].Vel[k];
        }
      local[n].State = BHP[b];
    }

  for(int task = 0; task < NTask; task++)
    {
      byte_counts[task] = counts[task] * (int)sizeof(struct bh_ffr_merge_summary);
      byte_displs[task] = displs[task] * (int)sizeof(struct bh_ffr_merge_summary);
    }

  MPI_Allgatherv(local, local_n * (int)sizeof(struct bh_ffr_merge_summary), MPI_BYTE,
                 all, byte_counts, byte_displs, MPI_BYTE, MPI_COMM_WORLD);

  qsort(all, nall, sizeof(*all), bh_ffr_merge_compare_id);

  int *parent = malloc(nall * sizeof(int));
  if(parent == NULL)
    terminate("BH_FFR: failed to allocate merger union-find state");
  for(int i = 0; i < nall; i++)
    parent[i] = i;

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  const double rmerge = 2.0 * All.BHAccretionRadius / a;
  const double rmerge2 = rmerge * rmerge;

  for(int i = 0; i < nall; i++)
    for(int j = i + 1; j < nall; j++)
      {
        MyDouble xtmp, ytmp, ztmp;
        const double dx = NEAREST_X(all[j].Pos[0] - all[i].Pos[0]);
        const double dy = NEAREST_Y(all[j].Pos[1] - all[i].Pos[1]);
        const double dz = NEAREST_Z(all[j].Pos[2] - all[i].Pos[2]);
        const double r2 = dx * dx + dy * dy + dz * dz;

        if(r2 <= rmerge2)
          bh_ffr_merge_union(parent, i, j);
      }

  for(int i = 0; i < nall; i++)
    parent[i] = bh_ffr_merge_find(parent, i);

  int global_mergers = 0;

  for(int root = 0; root < nall; root++)
    {
      if(parent[root] != root)
        continue;

      int nmembers = 0;
      for(int i = 0; i < nall; i++)
        if(parent[i] == root)
          nmembers++;

      if(nmembers < 2)
        continue;

      global_mergers += nmembers - 1;

      struct bh_ffr_particle_data merged;
      double merged_vel[3] = {0, 0, 0};
      int checked_members = 0;
      bh_ffr_build_merged_state(all, parent, root, nall, &merged, merged_vel, &checked_members);
      if(checked_members != nmembers)
        terminate("BH_FFR: merger component count mismatch %d != %d", checked_members, nmembers);

      const MyIDType survivor_id = all[root].ID;
      const int ps = bh_ffr_local_particle_from_id(survivor_id);

      if(ps >= 0)
        {
          const int bs = P[ps].BHDataIndex;
          double old_total_mass = 0.0;
          for(int i = 0; i < nall; i++)
            if(parent[i] == root)
              old_total_mass += all[i].DynMass;

          BHP[bs] = merged;
          P[ps].Mass = old_total_mass;
          for(int k = 0; k < 3; k++)
            P[ps].Vel[k] = merged_vel[k];

          printf("BH_FFR: merger survivor ID=%llu task=%d members=%d Mdyn=%g MBH=%g Mres=%g Mwindbuf=%g "
                 "Vx=%g Vy=%g Vz=%g\n",
                 (unsigned long long)survivor_id, ThisTask, nmembers, P[ps].Mass, BHP[bs].BHMass,
                 BHP[bs].ReservoirMass, BHP[bs].WindMassBuffer, P[ps].Vel[0], P[ps].Vel[1], P[ps].Vel[2]);
          fflush(stdout);
        }

      for(int i = 0; i < nall; i++)
        if(parent[i] == root && all[i].ID != survivor_id)
          {
            const int pl = bh_ffr_local_particle_from_id(all[i].ID);
            if(pl >= 0)
              {
                const MyIDType old_id = P[pl].ID;
                const double dx0 = nearest_x(all[i].Pos[0] - all[root].Pos[0]);
                const double dy0 = nearest_y(all[i].Pos[1] - all[root].Pos[1]);
                const double dz0 = nearest_z(all[i].Pos[2] - all[root].Pos[2]);
                const double separation = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0) * a;

                P[pl].Mass = 0.0;
                P[pl].ID = 0;
                P[pl].BHDataIndex = -1;
                P[pl].Type = BH_FFR_DM_PARTICLE_TYPE;
                for(int k = 0; k < 3; k++)
                  P[pl].Vel[k] = 0.0;

                printf("BH_FFR: merger consumed ID=%llu into ID=%llu task=%d separationProper=%g thresholdProper=%g\n",
                       (unsigned long long)old_id, (unsigned long long)survivor_id, ThisTask, separation,
                       2.0 * All.BHAccretionRadius);
                fflush(stdout);
              }
          }
    }

  if(global_mergers > 0)
    {
      timebin_cleanup_list_of_active_particles(&TimeBinsGravity);
      bh_ffr_rebuild_state_after_particle_changes();
      bh_ffr_validate_state("post-BH-merger");
    }

  free(parent);
  free(all);
  free(local);
  free(byte_displs);
  free(byte_counts);
  free(displs);
  free(counts);
}
