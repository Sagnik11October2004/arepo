#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Iteration 3 distributed gas-aperture search.
 *
 * BHAccretionRadius and BHFeedbackRadius are fixed proper lengths expressed
 * in AREPO's internal length unit. Particle coordinates are comoving when
 * ComovingIntegrationOn is enabled, so the tree-search radii are divided by
 * the current scale factor before entering the neighbour tree.
 *
 * This pass is intentionally read-only. It measures the gas environment and
 * supplies the nearby-gas timestep constraint. Gas removal, overlapping-BH
 * partitioning, angular-momentum accumulation, and feedback injection belong
 * to later iterations.
 */

static struct bh_ffr_gas_search_result *GasSearchResults;
static int GasSearchNTargets;

typedef struct
{
  MyDouble Pos[3];
  MyFloat SearchRadius;
  MyFloat AccretionRadius;
  MyFloat FeedbackRadius;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble AccretionMass;
  MyDouble FeedbackMass;
  long long AccretionCount;
  long long FeedbackCount;
  int MinHydroTimeBin;
} data_out;

static data_out *DataResult, *DataOut;

static double bh_ffr_proper_radius_to_search_radius(double proper_radius)
{
  if(!isfinite(proper_radius) || proper_radius <= 0)
    terminate("BH_FFR: invalid proper aperture radius %g", proper_radius);

  double fac = 1.0;
  if(All.ComovingIntegrationOn)
    {
      if(!isfinite(All.cf_atime) || All.cf_atime <= 0)
        terminate("BH_FFR: invalid scale factor cf_atime=%g while converting a proper aperture", All.cf_atime);
      fac = All.cf_atime;
    }

  const double radius = proper_radius / fac;
  if(!isfinite(radius) || radius <= 0)
    terminate("BH_FFR: invalid coordinate-space aperture radius %g from proper radius %g", radius, proper_radius);

  return radius;
}

static int bh_ffr_search_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= GasSearchNTargets)
    terminate("BH_FFR: gas-search target %d outside [0,%d) in %s", target, GasSearchNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active BH particle index %d for gas-search target %d in %s", p, target, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact state for gas-search particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: gas-search particle ID=%llu is not drifted to Ti_Current=%lld in %s", (unsigned long long)P[p].ID,
              (long long)All.Ti_Current, where);

  return p;
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_search_particle_from_target(target, "particle2in");

  in->Pos[0] = P[p].Pos[0];
  in->Pos[1] = P[p].Pos[1];
  in->Pos[2] = P[p].Pos[2];

  const double racc = bh_ffr_proper_radius_to_search_radius(All.BHAccretionRadius);
  const double rfb = bh_ffr_proper_radius_to_search_radius(All.BHFeedbackRadius);

  in->AccretionRadius = racc;
  in->FeedbackRadius = rfb;
  in->SearchRadius = dmax(racc, rfb);
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= GasSearchNTargets)
    terminate("BH_FFR: gas-search result target %d outside [0,%d)", target, GasSearchNTargets);

  struct bh_ffr_gas_search_result *res = &GasSearchResults[target];

  if(mode == MODE_LOCAL_PARTICLES)
    {
      res->AccretionMass = out->AccretionMass;
      res->FeedbackMass = out->FeedbackMass;
      res->AccretionCount = out->AccretionCount;
      res->FeedbackCount = out->FeedbackCount;
      res->MinHydroTimeBin = out->MinHydroTimeBin;
    }
  else
    {
      res->AccretionMass += out->AccretionMass;
      res->FeedbackMass += out->FeedbackMass;
      res->AccretionCount += out->AccretionCount;
      res->FeedbackCount += out->FeedbackCount;
      if(out->MinHydroTimeBin < res->MinHydroTimeBin)
        res->MinHydroTimeBin = out->MinHydroTimeBin;
    }
}

#include "../utils/generic_comm_helpers2.h"

static int bh_ffr_gas_search_evaluate(int target, int mode, int threadid);

static void kernel_local(void)
{
  const int threadid = get_thread_num();

  for(int task = 0; task < NTask; task++)
    Thread[threadid].Exportflag[task] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      const int target = NextParticle++;
      if(target >= GasSearchNTargets)
        break;

      bh_ffr_gas_search_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_gas_search_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static int bh_ffr_gas_search_evaluate(int target, int mode, int threadid)
{
  data_in local, *target_data;
  int numnodes, *firstnode;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  data_out out;
  memset(&out, 0, sizeof(out));
  out.MinHydroTimeBin = TIMEBINS;

  const double racc2 = target_data->AccretionRadius * target_data->AccretionRadius;
  const double rfb2 = target_data->FeedbackRadius * target_data->FeedbackRadius;

  const int nfound = ngb_treefind_variable_threads(target_data->Pos, target_data->SearchRadius, target, mode, threadid, numnodes,
                                                    firstnode);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_FFR: gas neighbour index %d outside NumGas=%d", j, NumGas);

      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      const double r2 = Thread[threadid].R2list[n];

      if(r2 <= racc2)
        {
          out.AccretionCount++;
          out.AccretionMass += P[j].Mass;

          int bin = P[j].TimeBinHydro;
          if(bin < 0)
            bin = -bin - 1;
          if(bin > 0 && bin < out.MinHydroTimeBin)
            out.MinHydroTimeBin = bin;
        }

      if(r2 <= rfb2)
        {
          out.FeedbackCount++;
          out.FeedbackMass += P[j].Mass;
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

void bh_ffr_collect_gas_environment(struct bh_ffr_gas_search_result *results)
{
  if(NumActiveBHFFR <= 0)
    return;
  if(results == NULL)
    terminate("BH_FFR: NULL result buffer for %d active BH gas searches", NumActiveBHFFR);

  GasSearchResults = results;
  GasSearchNTargets = NumActiveBHFFR;

  for(int n = 0; n < GasSearchNTargets; n++)
    {
      memset(&results[n], 0, sizeof(results[n]));
      results[n].MinHydroTimeBin = TIMEBINS;
      bh_ffr_search_particle_from_target(n, "bh_ffr_collect_gas_environment");
    }

  if(All.TotNumGas > 0)
    {
      generic_set_MaxNexport();
      generic_comm_pattern(GasSearchNTargets, kernel_local, kernel_imported);
    }

  for(int n = 0; n < GasSearchNTargets; n++)
    {
      if(results[n].MinHydroTimeBin == TIMEBINS)
        results[n].MinHydroTimeBin = -1;

      if(results[n].AccretionCount < 0 || results[n].FeedbackCount < 0 || !isfinite(results[n].AccretionMass) ||
         results[n].AccretionMass < 0 || !isfinite(results[n].FeedbackMass) || results[n].FeedbackMass < 0)
        terminate("BH_FFR: invalid distributed gas-search result for active target %d", n);
    }

  GasSearchResults = NULL;
  GasSearchNTargets = 0;
}

void bh_ffr_refresh_gas_neighbour_cache(void)
{
  if(NumActiveBHFFR <= 0)
    return;

  struct bh_ffr_gas_search_result *results =
      (struct bh_ffr_gas_search_result *)mymalloc("BHFFRGasSearchResults", NumActiveBHFFR * sizeof(*results));

  bh_ffr_collect_gas_environment(results);

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      bh_ffr_set_min_neighbour_timebin(p, results[n].MinHydroTimeBin);
    }

  myfree(results);
}

/*
 * Full-tree hook reserved for the cached local DM-dispersion query.
 *
 * The project specification requires a gravity-tree (not gas-tree) query and
 * BHDMNeighbours controls its eventual nearest-DM sample. The current project
 * has not defined a numeric particle type that universally means dark matter,
 * so Iteration 3 establishes and validates the full-tree execution point but
 * deliberately leaves BHP[].SigmaDM unchanged. A later iteration can add the
 * configured DM selector and velocity-moment calculation here without moving
 * the call site or querying an invalid gas-only tree.
 */
void bh_ffr_prepare_dm_environment_search(void)
{
  if(NumActiveBHFFR <= 0)
    return;

  if(Tree_MaxPart <= 0 || Tree_MaxNodes <= 0 || Nodes == NULL || Nextnode == NULL || Tree_Pos_list == NULL)
    terminate("BH_FFR: DM-environment hook called without a valid full gravity tree");

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
        terminate("BH_FFR: invalid active BH index %d in full-gravity-tree hook", p);
      if(P[p].Ti_Current != All.Ti_Current)
        terminate("BH_FFR: active BH particle ID=%llu is not drifted at full-gravity-tree hook", (unsigned long long)P[p].ID);
    }
}
