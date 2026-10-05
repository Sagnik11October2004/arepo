#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../fof/fof.h"
#include "../main/proto.h"

#if !defined(FOF)
#error "BLACKHOLE_FFR FoF seeding requires FOF"
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#error "BLACKHOLE_FFR FoF seeding requires FOF_SECONDARY_LINK_TYPES"
#endif

#if ((FOF_PRIMARY_LINK_TYPES & 2) == 0)
#error "BLACKHOLE_FFR FoF seeding expects Type-1 dark matter in FOF_PRIMARY_LINK_TYPES"
#endif

#if ((FOF_SECONDARY_LINK_TYPES & 1) == 0)
#error "BLACKHOLE_FFR FoF seeding requires gas (Type 0) in FOF_SECONDARY_LINK_TYPES"
#endif

#if ((FOF_SECONDARY_LINK_TYPES & 32) == 0)
#error "BLACKHOLE_FFR FoF seeding requires Type-5 BHs in FOF_SECONDARY_LINK_TYPES"
#endif

struct bh_ffr_seed_candidate
{
  int GrNr;
  double HaloMass;
  MyIDType SeedID;
  MyDouble SeedAxis[3];
  int HasCoherentAxis;
};

struct bh_ffr_seed_central
{
  MyIDType ID;
  MyDouble Pos[3];
  MyFloat Vel[3];
  double Density;
  double Mass;
  int Task;
  int Index;
};

static int SeedSelfTestDone;

static int bh_ffr_seed_candidate_compare(const void *a, const void *b)
{
  const struct bh_ffr_seed_candidate *ca = (const struct bh_ffr_seed_candidate *)a;
  const struct bh_ffr_seed_candidate *cb = (const struct bh_ffr_seed_candidate *)b;
  if(ca->GrNr < cb->GrNr)
    return -1;
  if(ca->GrNr > cb->GrNr)
    return 1;
  return 0;
}

static double bh_ffr_seed_msun_to_code_mass(double mass_msun)
{
  if(!isfinite(mass_msun) || !(mass_msun > 0) || !isfinite(All.UnitMass_in_g) || !(All.UnitMass_in_g > 0) ||
     !isfinite(All.HubbleParam) || !(All.HubbleParam > 0))
    terminate("BH_FFR: invalid seed mass conversion M=%g UnitMass=%g h=%g", mass_msun, All.UnitMass_in_g, All.HubbleParam);

  return mass_msun * SOLAR_MASS * All.HubbleParam / All.UnitMass_in_g;
}

static double bh_ffr_seed_distance2(const MyDouble pos[3], int i)
{
  const double dx = fof_periodic(P[i].Pos[0] - pos[0]);
  const double dy = fof_periodic(P[i].Pos[1] - pos[1]);
  const double dz = fof_periodic(P[i].Pos[2] - pos[2]);
  return dx * dx + dy * dy + dz * dz;
}

static void bh_ffr_seed_self_test(void)
{
  const double central = 2.0, target = 10.0, fmax = 0.5;
  const double donor_mass[3] = {5.0, 7.0, 9.0};
  double available = 0;
  for(int i = 0; i < 3; i++)
    available += fmax * donor_mass[i];

  const double need = target - central;
  if(!(need > 0) || available < need)
    terminate("BH_FFR: seed self-test failed availability");

  const double q = need / available;
  double removed = 0;
  for(int i = 0; i < 3; i++)
    removed += q * fmax * donor_mass[i];

  if(fabs((central + removed) - target) > 1.0e-14 * target)
    terminate("BH_FFR: seed self-test failed conservative target mass");
}

static struct bh_ffr_seed_candidate *bh_ffr_collect_seed_candidates(int *ncandidates, double halo_threshold_code)
{
  int nlocal = 0;
  for(int i = 0; i < Ngroups; i++)
    if(Group[i].Mass >= halo_threshold_code && Group[i].LenType[BH_FFR_PARTICLE_TYPE] == 0 && Group[i].LenType[0] > 0)
      nlocal++;

  struct bh_ffr_seed_candidate *local =
      (struct bh_ffr_seed_candidate *)mymalloc("BHFFRSeedCandidatesLocal", (nlocal > 0 ? nlocal : 1) * sizeof(*local));

  int n = 0;
  for(int i = 0; i < Ngroups; i++)
    if(Group[i].Mass >= halo_threshold_code && Group[i].LenType[BH_FFR_PARTICLE_TYPE] == 0 && Group[i].LenType[0] > 0)
      {
        local[n].GrNr = Group[i].GrNr;
        local[n].HaloMass = Group[i].Mass;
        local[n].SeedID = 0;
        local[n].HasCoherentAxis = 0;
        for(int k = 0; k < 3; k++)
          local[n].SeedAxis[k] = 0;
        n++;
      }

  int *counts = (int *)mymalloc("BHFFRSeedCandidateCounts", NTask * sizeof(int));
  int *displs = (int *)mymalloc("BHFFRSeedCandidateDispls", NTask * sizeof(int));
  int *counts_bytes = (int *)mymalloc("BHFFRSeedCandidateCountsBytes", NTask * sizeof(int));
  int *displs_bytes = (int *)mymalloc("BHFFRSeedCandidateDisplsBytes", NTask * sizeof(int));

  MPI_Allgather(&nlocal, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  int ntotal = 0;
  for(int task = 0; task < NTask; task++)
    {
      if(counts[task] > INT_MAX / (int)sizeof(*local))
        terminate("BH_FFR: too many FoF seed candidates on task %d", task);
      displs[task] = ntotal;
      ntotal += counts[task];
      counts_bytes[task] = counts[task] * (int)sizeof(*local);
      displs_bytes[task] = displs[task] * (int)sizeof(*local);
    }

  /* The returned candidate list outlives the scratch arrays below. Allocate
   * it outside AREPO's LIFO arena so the scratch blocks can be released here. */
  struct bh_ffr_seed_candidate *all =
      (struct bh_ffr_seed_candidate *)malloc((ntotal > 0 ? ntotal : 1) * sizeof(*all));
  if(all == NULL)
    terminate("BH_FFR: failed to allocate %d global FoF seed candidates", ntotal);

  MPI_Allgatherv(local, nlocal * (int)sizeof(*local), MPI_BYTE, all, counts_bytes, displs_bytes, MPI_BYTE, MPI_COMM_WORLD);

  if(ntotal > 1)
    qsort(all, ntotal, sizeof(*all), bh_ffr_seed_candidate_compare);

  for(int i = 1; i < ntotal; i++)
    if(all[i].GrNr == all[i - 1].GrNr)
      terminate("BH_FFR: duplicate global FoF group number %d in seed candidate list", all[i].GrNr);

  myfree(displs_bytes);
  myfree(counts_bytes);
  myfree(displs);
  myfree(counts);
  myfree(local);

  *ncandidates = ntotal;
  return all;
}

static int bh_ffr_find_seed_central(int grnr, double seed_mass_code, struct bh_ffr_seed_central *central)
{
  double local_density = -DBL_MAX;
  unsigned long long local_id = ULLONG_MAX;
  int local_index = -1;

  for(int i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0 && PS[i].GrNr == grnr)
      {
        const double density = SphP[i].Density;
        const unsigned long long id = (unsigned long long)P[i].ID;
        if(isfinite(density) && (density > local_density || (density == local_density && id < local_id)))
          {
            local_density = density;
            local_id = id;
            local_index = i;
          }
      }

  double global_density = -DBL_MAX;
  MPI_Allreduce(&local_density, &global_density, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(global_density == -DBL_MAX)
    return 0;

  unsigned long long tie_id = (local_density == global_density) ? local_id : ULLONG_MAX;
  unsigned long long global_id = ULLONG_MAX;
  MPI_Allreduce(&tie_id, &global_id, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);

  int owner_local = (local_index >= 0 && (unsigned long long)P[local_index].ID == global_id &&
                     SphP[local_index].Density == global_density) ? ThisTask : INT_MAX;
  int owner = INT_MAX;
  MPI_Allreduce(&owner_local, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if(owner == INT_MAX)
    terminate("BH_FFR: failed to locate owner of FoF seed central ID=%llu", global_id);

  memset(central, 0, sizeof(*central));
  if(ThisTask == owner)
    {
      central->ID = P[local_index].ID;
      central->Density = SphP[local_index].Density;
      central->Mass = P[local_index].Mass;
      central->Task = ThisTask;
      central->Index = local_index;
      for(int k = 0; k < 3; k++)
        {
          central->Pos[k] = P[local_index].Pos[k];
          central->Vel[k] = P[local_index].Vel[k];
        }
    }

  MPI_Bcast(central, (int)sizeof(*central), MPI_BYTE, owner, MPI_COMM_WORLD);

  if(!(central->Mass > 0) || !isfinite(central->Density))
    terminate("BH_FFR: invalid FoF seed central ID=%llu mass=%g density=%g", (unsigned long long)central->ID,
              central->Mass, central->Density);

  /* The seed mass must exceed the mass of the actual densest gas cell.
   * Do not silently move the seed to a less-dense cell if this is violated. */
  if(!(central->Mass < seed_mass_code))
    return 0;

  return 1;
}

static int *bh_ffr_build_local_donor_list(int grnr, MyIDType central_id, int *ndonors)
{
  int n = 0;
  for(int i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].ID != central_id && P[i].Mass > 0 && PS[i].GrNr == grnr)
      n++;

  int *list = (int *)mymalloc("BHFFRSeedDonorList", (n > 0 ? n : 1) * sizeof(int));
  int j = 0;
  for(int i = 0; i < NumGas; i++)
    if(P[i].Type == 0 && P[i].ID != 0 && P[i].ID != central_id && P[i].Mass > 0 && PS[i].GrNr == grnr)
      list[j++] = i;

  *ndonors = n;
  return list;
}

static double bh_ffr_seed_measure_local_rotation(const struct bh_ffr_seed_central *central,
                                                   const int *donors, int ndonors, MyDouble axis[3])
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  const double rcoord = All.BHAccretionRadius / a;
  const double r2max = rcoord * rcoord;
  double local_mass = 0.0, local_c[3] = {0, 0, 0};

  for(int n = 0; n < ndonors; n++)
    {
      const int i = donors[n];
      if(bh_ffr_seed_distance2(central->Pos, i) > r2max)
        continue;
      double xtmp, ytmp, ztmp;
      const double dr[3] = {NEAREST_X(P[i].Pos[0] - central->Pos[0]),
                            NEAREST_Y(P[i].Pos[1] - central->Pos[1]),
                            NEAREST_Z(P[i].Pos[2] - central->Pos[2])};
      const double dv[3] = {P[i].Vel[0] - central->Vel[0], P[i].Vel[1] - central->Vel[1], P[i].Vel[2] - central->Vel[2]};
      const double ell[3] = {dr[1]*dv[2]-dr[2]*dv[1], dr[2]*dv[0]-dr[0]*dv[2], dr[0]*dv[1]-dr[1]*dv[0]};
      const double en = sqrt(ell[0]*ell[0] + ell[1]*ell[1] + ell[2]*ell[2]);
      if(!(en > 0) || !isfinite(en))
        continue;
      local_mass += P[i].Mass;
      for(int k = 0; k < 3; k++)
        local_c[k] += P[i].Mass * ell[k] / en;
    }

  double mass = 0.0, c[3] = {0, 0, 0};
  MPI_Allreduce(&local_mass, &mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_c, c, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double cn = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  const double coherence = mass > 0 ? cn / mass : 0.0;
  if(!isfinite(coherence) || coherence < 0 || coherence > 1.0 + 1.0e-10)
    terminate("BH_FFR: invalid local seed coherence=%g", coherence);
  for(int k = 0; k < 3; k++)
    axis[k] = (cn > 0 && coherence > All.BHMinCoherence) ? c[k] / cn : 0.0;
  return coherence;
}

static int bh_ffr_seed_one_candidate(struct bh_ffr_seed_candidate *candidate, double seed_mass_code)
{
  struct bh_ffr_seed_central central;
  if(!bh_ffr_find_seed_central(candidate->GrNr, seed_mass_code, &central))
    return 0;

  const double need = seed_mass_code - central.Mass;
  if(!(need > 0))
    return 0;

  int ndonors = 0;
  int *donors = bh_ffr_build_local_donor_list(candidate->GrNr, central.ID, &ndonors);

  double local_available = 0, local_r2max = 0;
  for(int n = 0; n < ndonors; n++)
    {
      const int i = donors[n];
      local_available += All.BHSeedMaxDonorFraction * P[i].Mass;
      local_r2max = dmax(local_r2max, bh_ffr_seed_distance2(central.Pos, i));
    }

  double available = 0, r2max = 0;
  MPI_Allreduce(&local_available, &available, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_r2max, &r2max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  const double mass_tol = 1.0e-12 * dmax(fabs(seed_mass_code), 1.0e-30);
  if(available + mass_tol < need)
    {
      myfree(donors);
      return 0;
    }

  double r2lo = 0, r2hi = r2max;
  if(r2hi > 0)
    for(int iter = 0; iter < 40; iter++)
      {
        const double r2mid = 0.5 * (r2lo + r2hi);
        double local_inside = 0;
        for(int n = 0; n < ndonors; n++)
          {
            const int i = donors[n];
            if(bh_ffr_seed_distance2(central.Pos, i) <= r2mid)
              local_inside += All.BHSeedMaxDonorFraction * P[i].Mass;
          }
        double inside_mid = 0;
        MPI_Allreduce(&local_inside, &inside_mid, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(inside_mid >= need)
          r2hi = r2mid;
        else
          r2lo = r2mid;
      }

  if(r2hi > 0)
    r2hi = nextafter(r2hi, DBL_MAX);

  double local_inside = 0;
  for(int n = 0; n < ndonors; n++)
    {
      const int i = donors[n];
      if(bh_ffr_seed_distance2(central.Pos, i) <= r2hi)
        local_inside += All.BHSeedMaxDonorFraction * P[i].Mass;
    }

  double inside = 0;
  MPI_Allreduce(&local_inside, &inside, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(inside + mass_tol < need)
    {
      r2hi = r2max;
      inside = available;
    }

  if(!(inside > 0) || inside + mass_tol < need)
    terminate("BH_FFR: invalid donor mass inside seed sphere group=%d inside=%g need=%g", candidate->GrNr, inside, need);

  double q = need / inside;
  if(q > 1.0 && q < 1.0 + 1.0e-12)
    q = 1.0;
  if(!(q > 0) || q > 1.0 || !isfinite(q))
    terminate("BH_FFR: invalid seed donor normalization q=%g group=%d", q, candidate->GrNr);

  candidate->SeedID = central.ID;
  const double seed_coherence = bh_ffr_seed_measure_local_rotation(&central, donors, ndonors, candidate->SeedAxis);
  candidate->HasCoherentAxis = seed_coherence > All.BHMinCoherence &&
                               (candidate->SeedAxis[0] != 0 || candidate->SeedAxis[1] != 0 || candidate->SeedAxis[2] != 0);

  double local_removed_mass = 0;
  double local_removed_momentum[3] = {0, 0, 0};

  if(ThisTask == central.Task)
    {
      const int i = central.Index;
      if(i < 0 || i >= NumGas || P[i].Type != 0 || P[i].ID != central.ID)
        terminate("BH_FFR: central gas cell moved before seed commit ID=%llu", (unsigned long long)central.ID);
      local_removed_mass += P[i].Mass;
      for(int k = 0; k < 3; k++)
        local_removed_momentum[k] += SphP[i].Momentum[k];
    }

  for(int n = 0; n < ndonors; n++)
    {
      const int i = donors[n];
      if(bh_ffr_seed_distance2(central.Pos, i) > r2hi)
        continue;

      const double oldmass = P[i].Mass;
      const double dm = q * All.BHSeedMaxDonorFraction * oldmass;
      if(!(dm >= 0) || dm >= oldmass)
        terminate("BH_FFR: unsafe seed donor removal ID=%llu oldmass=%g dm=%g", (unsigned long long)P[i].ID, oldmass, dm);

      const double keep = (oldmass - dm) / oldmass;
      local_removed_mass += dm;
      for(int k = 0; k < 3; k++)
        {
          local_removed_momentum[k] += (1.0 - keep) * SphP[i].Momentum[k];
          SphP[i].Momentum[k] *= keep;
        }

      P[i].Mass *= keep;
      SphP[i].Energy *= keep;
#ifdef PASSIVE_SCALARS
      for(int k = 0; k < PASSIVE_SCALARS; k++)
        SphP[i].PConservedScalars[k] *= keep;
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
      SphP[i].HighResMass *= keep;
#endif
      if(!(SphP[i].Volume > 0))
        terminate("BH_FFR: non-positive donor volume for gas ID=%llu", (unsigned long long)P[i].ID);
      SphP[i].Density = P[i].Mass / SphP[i].Volume;
      SphP[i].OldMass = P[i].Mass;
      set_pressure_of_cell(i);
      if(!(P[i].Mass > 0) || !isfinite(P[i].Mass) || !(SphP[i].Density > 0) || !isfinite(SphP[i].Density))
        terminate("BH_FFR: invalid donor state after seed removal ID=%llu", (unsigned long long)P[i].ID);
    }

  myfree(donors);

  double removed_mass = 0;
  double removed_momentum[3] = {0, 0, 0};
  MPI_Allreduce(&local_removed_mass, &removed_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_removed_momentum, removed_momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(fabs(removed_mass - seed_mass_code) > 1.0e-9 * dmax(fabs(seed_mass_code), 1.0e-30))
    terminate("BH_FFR: seed mass ledger failed group=%d removed=%g target=%g", candidate->GrNr, removed_mass, seed_mass_code);

  if(ThisTask == central.Task)
    {
      const int i = central.Index;
      P[i].Type = BH_FFR_PARTICLE_TYPE;
      P[i].SofteningType = All.SofteningTypeOfPartType[BH_FFR_PARTICLE_TYPE];
      P[i].Mass = seed_mass_code;
      P[i].BHDataIndex = -1;
      for(int k = 0; k < 3; k++)
        P[i].Vel[k] = removed_momentum[k] / seed_mass_code;
      if(!isfinite(P[i].Vel[0]) || !isfinite(P[i].Vel[1]) || !isfinite(P[i].Vel[2]))
        terminate("BH_FFR: non-finite seed velocity for ID=%llu", (unsigned long long)P[i].ID);
    }

  mpi_printf("BH_FFR: seeded FoF group %d at z=%g with Mhalo=%g Msun, MBH=%g Msun, central gas ID=%llu.\n",
             candidate->GrNr, All.cf_redshift,
             candidate->HaloMass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS), All.BHSeedMassMsun,
             (unsigned long long)central.ID);
  mpi_printf("BH_FFR: seed axis ID=%llu axis_source=%s coherence=%g.\n",
             (unsigned long long)central.ID, candidate->HasCoherentAxis ? "local-gas" : "deterministic-fallback",
             seed_coherence);
  return 1;
}

static int bh_ffr_rearrange_seeded_gas_cells(void)
{
  int converted = 0;
  for(int i = 0; i < NumGas; i++)
    if(P[i].Type != 0)
      {
        if(P[i].Type != BH_FFR_PARTICLE_TYPE)
          terminate("BH_FFR: unexpected non-gas Type-%d particle inside gas block during seed rearrangement", P[i].Type);

        struct particle_data psave = P[i];
        const peanokey keysave = Key[i];
        P[i] = P[NumGas - 1];
        SphP[i] = SphP[NumGas - 1];
        Key[i] = Key[NumGas - 1];
        P[NumGas - 1] = psave;
        Key[NumGas - 1] = keysave;
        NumGas--;
        converted++;
        i--;
      }
  return converted;
}

int bh_ffr_seed_from_fof(void)
{
#ifdef MHD
  terminate("BH_FFR: FoF seed formation is not yet enabled with MHD; a magnetic-flux conversion policy is required first");
#endif

  if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
    terminate("BH_FFR: FoF seeding called away from a full synchronization point");
  if(!All.ComovingIntegrationOn)
    terminate("BH_FFR: redshift-threshold FoF seeding requires ComovingIntegrationOn");

  if(!SeedSelfTestDone)
    {
      bh_ffr_seed_self_test();
      SeedSelfTestDone = 1;
    }

  if(!isfinite(All.BHSeedHaloMassMsun) || !(All.BHSeedHaloMassMsun > 0) || !isfinite(All.BHSeedMassMsun) ||
     !(All.BHSeedMassMsun > 0) || !isfinite(All.BHSeedMinRedshift) || All.BHSeedMinRedshift < 0 ||
     !isfinite(All.BHSeedMaxDonorFraction) || !(All.BHSeedMaxDonorFraction > 0) || !(All.BHSeedMaxDonorFraction < 1))
    terminate("BH_FFR: invalid seed parameters halo=%g Msun seed=%g Msun zmin=%g donorfrac=%g", All.BHSeedHaloMassMsun,
              All.BHSeedMassMsun, All.BHSeedMinRedshift, All.BHSeedMaxDonorFraction);

  if(!(All.cf_redshift > All.BHSeedMinRedshift))
    return 0;

  const double halo_threshold_code = bh_ffr_seed_msun_to_code_mass(All.BHSeedHaloMassMsun);
  const double seed_mass_code = bh_ffr_seed_msun_to_code_mass(All.BHSeedMassMsun);

  int ncandidates = 0;
  struct bh_ffr_seed_candidate *candidates = bh_ffr_collect_seed_candidates(&ncandidates, halo_threshold_code);

  int seeded_transactions = 0;
  for(int c = 0; c < ncandidates; c++)
    seeded_transactions += bh_ffr_seed_one_candidate(&candidates[c], seed_mass_code);

  int local_converted = bh_ffr_rearrange_seeded_gas_cells();
  int global_converted = 0;
  MPI_Allreduce(&local_converted, &global_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(global_converted != seeded_transactions)
    terminate("BH_FFR: seed conversion count mismatch transactions=%d converted=%d", seeded_transactions, global_converted);
  if(global_converted == 0)
    {
      free(candidates);
      return 0;
    }

  All.TotNumGas -= global_converted;
  if(All.TotNumGas < 0)
    terminate("BH_FFR: negative global gas count after seeding");

  bh_ffr_rebuild_state_after_particle_changes();
  for(int c = 0; c < ncandidates; c++)
    if(candidates[c].SeedID != 0 && candidates[c].HasCoherentAxis)
      for(int b = 0; b < NumBHFFR; b++)
        if(BHP[b].ParticleID == candidates[c].SeedID)
          for(int k = 0; k < 3; k++)
            {
              BHP[b].DiscDir[k] = candidates[c].SeedAxis[k];
              BHP[b].JetDir[k] = candidates[c].SeedAxis[k];
            }
  free(candidates);
  reconstruct_timebins();

  /* The FoF caller still owns PS and other movable arena blocks here. The
   * gas neighbour tree is stale after the gas-to-BH conversion, but it is
   * rebuilt by fof_fof() only after FoF has released those scratch blocks. */
  bh_ffr_validate_state("FoF seeding");

  long long local_bh = NumBHFFR, global_bh = 0;
  MPI_Allreduce(&local_bh, &global_bh, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("BH_FFR: FoF full-step seeding created %d seed(s); total Type-%d BH count is now %lld.\n",
             global_converted, BH_FFR_PARTICLE_TYPE, global_bh);
  return global_converted;
}
