#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../gravity/forcetree.h"
#include "../main/proto.h"

/*
 * Iteration 10: cached local dark-matter velocity dispersion.
 *
 * The public neighbour tree is gas-only, so this module walks AREPO's full
 * gravity tree at synchronization points where that tree contains the complete
 * particle set.  It keeps the globally nearest BHDMNeighbours Type-1 dark
 * matter particles around every active BH and computes the one-dimensional
 * physical peculiar velocity dispersion about their local mean:
 *
 *   sigma_DM^2 = < |v - <v>|^2 > / 3 .
 *
 * SigmaDM is then cached in BHP[] until the next full-tree update.  No extra
 * persistent cache metadata is added, preserving the native BHP restart ABI.
 *
 * If fewer than BHDMNeighbours DM particles exist globally for a target, the
 * same gravity-tree query is repeated for gas and its nearest sample is used
 * as a fallback.  The fallback is explicitly logged; no persistent fallback
 * flag is added because the current BHP layout is restart-frozen.
 */

#define BH_FFR_DM_MAX_NEIGHBOURS 256
#define BH_FFR_GAS_PARTICLE_TYPE 0

/* Iteration-11 deliberately keeps these simple and compile-time. They can be
 * promoted to runtime parameters after physical convergence tests without
 * changing the persistent BHP restart layout. */
#define BH_FFR_DF_COULOMB_LOG 3.0
#define BH_FFR_DF_MAX_DAMPING_FRACTION 0.5

struct bh_ffr_dm_candidate
{
  MyDouble R2;
  MyDouble Mass;
  MyFloat Vel[3];
};

struct bh_ffr_dm_result
{
  int Count;
  struct bh_ffr_dm_candidate C[BH_FFR_DM_MAX_NEIGHBOURS];
};

static struct bh_ffr_dm_result *DMResults;
static int DMNTargets;
static int DMQueryType;

typedef struct
{
  MyDouble Pos[3];
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  int Count;
  struct bh_ffr_dm_candidate C[BH_FFR_DM_MAX_NEIGHBOURS];
} data_out;

static data_out *DataResult, *DataOut;

static int bh_ffr_dm_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= DMNTargets)
    terminate("BH_FFR: DM target=%d outside [0,%d) in %s", target, DMNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active BH particle=%d in DM query %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact DM-query state for ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: BH ID=%llu is not drifted at DM query in %s", (unsigned long long)P[p].ID, where);

  return p;
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_dm_particle_from_target(target, "particle2in");
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[p].Pos[k];
  in->Firstnode = firstnode;
}

static void bh_ffr_dm_insert(struct bh_ffr_dm_result *res, double r2, double mass, const MyFloat vel[3])
{
  const int want = All.BHDMNeighbours;
  if(want < 1 || want > BH_FFR_DM_MAX_NEIGHBOURS)
    terminate("BH_FFR: BHDMNeighbours=%d outside supported [1,%d]", want, BH_FFR_DM_MAX_NEIGHBOURS);

  if(!isfinite(r2) || r2 < 0 || !isfinite(mass) || mass < 0)
    terminate("BH_FFR: invalid DM-neighbour r2=%g mass=%g", r2, mass);

  int pos;

  if(res->Count < want)
    {
      pos = res->Count++;
    }
  else
    {
      if(r2 >= res->C[want - 1].R2)
        return;
      pos = want - 1;
    }

  res->C[pos].R2 = r2;
  res->C[pos].Mass = mass;
  for(int k = 0; k < 3; k++)
    {
      if(!isfinite(vel[k]))
        terminate("BH_FFR: non-finite DM-neighbour velocity component=%g", (double)vel[k]);
      res->C[pos].Vel[k] = vel[k];
    }

  while(pos > 0 && res->C[pos].R2 < res->C[pos - 1].R2)
    {
      struct bh_ffr_dm_candidate tmp = res->C[pos - 1];
      res->C[pos - 1] = res->C[pos];
      res->C[pos] = tmp;
      pos--;
    }
}

static void bh_ffr_dm_merge(struct bh_ffr_dm_result *dst, const data_out *src)
{
  if(src->Count < 0 || src->Count > All.BHDMNeighbours || src->Count > BH_FFR_DM_MAX_NEIGHBOURS)
    terminate("BH_FFR: invalid imported DM-neighbour count=%d", src->Count);

  for(int n = 0; n < src->Count; n++)
    bh_ffr_dm_insert(dst, src->C[n].R2, src->C[n].Mass, src->C[n].Vel);
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= DMNTargets)
    terminate("BH_FFR: DM result target=%d outside [0,%d)", target, DMNTargets);

  if(mode == MODE_LOCAL_PARTICLES)
    {
      DMResults[target].Count = 0;
      bh_ffr_dm_merge(&DMResults[target], out);
    }
  else
    bh_ffr_dm_merge(&DMResults[target], out);
}

#include "../utils/generic_comm_helpers2.h"

static int bh_ffr_dm_evaluate(int target, int mode, int threadid);

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
      if(target >= DMNTargets)
        break;

      bh_ffr_dm_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_dm_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static double bh_ffr_dm_node_min_r2(const struct NODE *node, const MyDouble pos[3])
{
  MyDouble xtmp, ytmp, ztmp;
  const double half = 0.5 * node->len;

  double dx = fabs(GRAVITY_NEAREST_X(node->center[0] - pos[0])) - half;
  double dy = fabs(GRAVITY_NEAREST_Y(node->center[1] - pos[1])) - half;
  double dz = fabs(GRAVITY_NEAREST_Z(node->center[2] - pos[2])) - half;

  if(dx < 0)
    dx = 0;
  if(dy < 0)
    dy = 0;
  if(dz < 0)
    dz = 0;

  return dx * dx + dy * dy + dz * dz;
}

static void bh_ffr_dm_consider(struct bh_ffr_dm_result *out, const MyDouble pos[3], const MyDouble particle_pos[3],
                               double mass, const MyFloat stored_vel[3])
{
  MyDouble xtmp, ytmp, ztmp;
  const double dx = GRAVITY_NEAREST_X(particle_pos[0] - pos[0]);
  const double dy = GRAVITY_NEAREST_Y(particle_pos[1] - pos[1]);
  const double dz = GRAVITY_NEAREST_Z(particle_pos[2] - pos[2]);
  const double r2 = dx * dx + dy * dy + dz * dz;

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_FFR: invalid scale factor=%g in DM velocity query", a);

  MyFloat vel[3];
  for(int k = 0; k < 3; k++)
    vel[k] = stored_vel[k] / a;

  bh_ffr_dm_insert(out, r2, mass, vel);
}

static int bh_ffr_dm_evaluate(int target, int mode, int threadid)
{
  data_in local, *in;
  int numnodes, *firstnode;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  struct bh_ffr_dm_result tmp;
  memset(&tmp, 0, sizeof(tmp));

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == MODE_LOCAL_PARTICLES)
        no = Tree_MaxPart;
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart)
            {
              const int p = no;
              no = Nextnode[no];

              if(P[p].Type != DMQueryType || P[p].ID == 0 || !(P[p].Mass > 0))
                continue;

              MyFloat vel[3] = {P[p].Vel[0], P[p].Vel[1], P[p].Vel[2]};
              bh_ffr_dm_consider(&tmp, in->Pos, &Tree_Pos_list[3 * p], P[p].Mass, vel);
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)
            {
              if(mode == MODE_IMPORTED_PARTICLES && no < Tree_FirstNonTopLevelNode)
                break;

              struct NODE *node = &Nodes[no];
              const int sibling = node->u.d.sibling;

              if(tmp.Count >= All.BHDMNeighbours)
                {
                  const double cutoff2 = tmp.C[All.BHDMNeighbours - 1].R2;
                  if(bh_ffr_dm_node_min_r2(node, in->Pos) > cutoff2)
                    {
                      no = sibling;
                      continue;
                    }
                }

              no = node->u.d.nextnode;
            }
          else if(no >= Tree_ImportedNodeOffset)
            {
              const int n = no - Tree_ImportedNodeOffset;
              no = Nextnode[no - Tree_MaxNodes];

              if(n < 0 || n >= Tree_NumPartImported)
                terminate("BH_FFR: imported gravity point=%d outside [0,%d)", n, Tree_NumPartImported);

              if(Tree_Points[n].Type != DMQueryType || !(Tree_Points[n].Mass > 0))
                continue;

              bh_ffr_dm_consider(&tmp, in->Pos, Tree_Points[n].Pos, Tree_Points[n].Mass, Tree_Points[n].Vel);
            }
          else
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("BH_FFR: pseudo particle encountered in imported DM tree walk");

              tree_treefind_export_node_threads(no, target, threadid);
              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  data_out out;
  memset(&out, 0, sizeof(out));
  out.Count = tmp.Count;
  for(int n = 0; n < tmp.Count; n++)
    out.C[n] = tmp.C[n];

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static void bh_ffr_dm_run_query(int particle_type, struct bh_ffr_dm_result *results)
{
  DMQueryType = particle_type;
  DMResults = results;
  DMNTargets = NumActiveBHFFR;

  for(int n = 0; n < DMNTargets; n++)
    {
      memset(&results[n], 0, sizeof(results[n]));
      bh_ffr_dm_particle_from_target(n, "bh_ffr_dm_run_query");
    }

  generic_set_MaxNexport();
  generic_comm_pattern(DMNTargets, kernel_local, kernel_imported);

  DMResults = NULL;
  DMNTargets = 0;
}

static double bh_ffr_dm_sigma_1d(const struct bh_ffr_dm_result *res)
{
  if(res->Count < 2)
    return 0.0;

  double mean[3] = {0, 0, 0};
  for(int n = 0; n < res->Count; n++)
    for(int k = 0; k < 3; k++)
      mean[k] += res->C[n].Vel[k];

  for(int k = 0; k < 3; k++)
    mean[k] /= res->Count;

  double sum2 = 0.0;
  for(int n = 0; n < res->Count; n++)
    for(int k = 0; k < 3; k++)
      {
        const double dv = res->C[n].Vel[k] - mean[k];
        sum2 += dv * dv;
      }

  double sigma2 = sum2 / (3.0 * res->Count);
  if(sigma2 < 0 && sigma2 > -1.0e-14 * fmax(1.0, sum2))
    sigma2 = 0;

  if(!isfinite(sigma2) || sigma2 < 0)
    terminate("BH_FFR: invalid local velocity variance=%g", sigma2);

  return sqrt(sigma2);
}

static void bh_ffr_dm_mean_velocity(const struct bh_ffr_dm_result *res, double mean[3])
{
  for(int k = 0; k < 3; k++)
    mean[k] = 0.0;

  if(res->Count <= 0)
    return;

  for(int n = 0; n < res->Count; n++)
    for(int k = 0; k < 3; k++)
      mean[k] += res->C[n].Vel[k];

  for(int k = 0; k < 3; k++)
    mean[k] /= res->Count;
}

static double bh_ffr_dm_density_code(const struct bh_ffr_dm_result *res)
{
  if(res->Count < 2)
    return 0.0;

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  const double r = sqrt(res->C[res->Count - 1].R2) * a;
  if(!isfinite(r) || !(r > 0))
    return 0.0;

  double mass = 0.0;
  for(int n = 0; n < res->Count; n++)
    mass += res->C[n].Mass;

  const double volume = (4.0 * M_PI / 3.0) * r * r * r;
  const double rho = mass / volume;

  if(!isfinite(rho) || rho < 0)
    terminate("BH_FFR: invalid local DM density=%g from mass=%g r=%g", rho, mass, r);

  return rho;
}

static void bh_ffr_apply_dynamical_friction(int p, const struct bh_ffr_dm_result *dm)
{
  if(dm->Count < 2 || P[p].Mass <= 0)
    return;

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid BH state in dynamical friction for ID=%llu", (unsigned long long)P[p].ID);

  const double dt = bh_ffr_get_elapsed_time_code_time(p);
  if(dt <= 0)
    return;

  const double sigma = bh_ffr_dm_sigma_1d(dm);
  const double rho = bh_ffr_dm_density_code(dm);
  if(!(sigma > 0) || !(rho > 0))
    return;

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  double mean[3], vrel[3];
  bh_ffr_dm_mean_velocity(dm, mean);

  double v2 = 0.0;
  for(int k = 0; k < 3; k++)
    {
      vrel[k] = P[p].Vel[k] / a - mean[k];
      v2 += vrel[k] * vrel[k];
    }

  if(!(v2 > 0) || !isfinite(v2))
    return;

  const double v = sqrt(v2);
  const double x = v / (M_SQRT2 * sigma);
  double fx;

  if(x < 1.0e-3)
    fx = 4.0 * x * x * x / (3.0 * sqrt(M_PI));
  else
    fx = erf(x) - 2.0 * x * exp(-x * x) / sqrt(M_PI);

  if(fx < 0 && fx > -1.0e-14)
    fx = 0.0;
  if(!isfinite(fx) || fx < 0)
    terminate("BH_FFR: invalid Chandrasekhar F(X)=%g for X=%g", fx, x);
  if(fx == 0)
    return;

  const double adf = 4.0 * M_PI * All.G * All.G * P[p].Mass * rho * BH_FFR_DF_COULOMB_LOG * fx / v2;
  if(!isfinite(adf) || adf < 0)
    terminate("BH_FFR: invalid dynamical-friction acceleration=%g", adf);
  if(adf == 0)
    return;

  const double tdf = v / adf;
  double frac = -expm1(-dt / tdf);
  if(frac > BH_FFR_DF_MAX_DAMPING_FRACTION)
    frac = BH_FFR_DF_MAX_DAMPING_FRACTION;

  if(!isfinite(frac) || frac < 0 || frac > BH_FFR_DF_MAX_DAMPING_FRACTION)
    terminate("BH_FFR: invalid dynamical-friction damping fraction=%g", frac);

  const double v_before = v;
  for(int k = 0; k < 3; k++)
    P[p].Vel[k] -= a * frac * vrel[k];

  double vafter2 = 0.0;
  for(int k = 0; k < 3; k++)
    {
      const double dv = P[p].Vel[k] / a - mean[k];
      vafter2 += dv * dv;
    }
  const double v_after = sqrt(vafter2);

  if(v_after > v_before * (1.0 + 2.0e-12))
    terminate("BH_FFR: dynamical friction increased relative velocity for ID=%llu", (unsigned long long)P[p].ID);

  printf("BH_FFR: dynamical friction ID=%llu task=%d rhoDM=%g sigmaDM=%g vrel0=%g vrel1=%g "
         "tdfCode=%g dtCode=%g frac=%g lnLambda=%g\n",
         (unsigned long long)P[p].ID, ThisTask, rho, sigma, v_before, v_after, tdf, dt, frac,
         BH_FFR_DF_COULOMB_LOG);
  fflush(stdout);
}

static void bh_ffr_dm_self_test(void)
{
  struct bh_ffr_dm_result res;
  memset(&res, 0, sizeof(res));

  const MyFloat v0[3] = {-1, 0, 0};
  const MyFloat v1[3] = {1, 0, 0};
  const MyFloat v2[3] = {0, 0, 0};

  bh_ffr_dm_insert(&res, 4.0, 1.0, v0);
  bh_ffr_dm_insert(&res, 1.0, 1.0, v1);
  bh_ffr_dm_insert(&res, 2.0, 1.0, v2);

  if(res.Count != 3 || res.C[0].R2 != 1.0 || res.C[1].R2 != 2.0 || res.C[2].R2 != 4.0)
    terminate("BH_FFR: DM nearest-neighbour ordering self-test failed");

  const double sigma = bh_ffr_dm_sigma_1d(&res);
  const double expected = sqrt(2.0 / 9.0);
  if(fabs(sigma - expected) > 1.0e-13)
    terminate("BH_FFR: DM dispersion self-test failed sigma=%g expected=%g", sigma, expected);
}

void bh_ffr_prepare_dm_environment_search(void)
{
#ifdef HIERARCHICAL_GRAVITY
  /* The hook can also be reached with a partial hierarchical tree. SigmaDM is
   * refreshed only when the highest occupied gravity bin is synchronized,
   * i.e. when the complete gravity tree is present. */
  if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
    return;
#endif

  bh_ffr_build_active_list();

  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

  if(Tree_MaxPart <= 0 || Tree_MaxNodes <= 0 || Nodes == NULL || Nextnode == NULL || Tree_Pos_list == NULL)
    terminate("BH_FFR: SigmaDM query called without a valid full gravity tree");

  if(All.BHDMNeighbours < 2 || All.BHDMNeighbours > BH_FFR_DM_MAX_NEIGHBOURS)
    terminate("BH_FFR: BHDMNeighbours=%d must lie in [2,%d]", All.BHDMNeighbours, BH_FFR_DM_MAX_NEIGHBOURS);

  static int self_test_done = 0;
  if(!self_test_done)
    {
      bh_ffr_dm_self_test();
      self_test_done = 1;
    }

  struct bh_ffr_dm_result *dm = (struct bh_ffr_dm_result *)mymalloc(
      "BHFFRDMResults", (NumActiveBHFFR > 0 ? NumActiveBHFFR : 1) * sizeof(*dm));

  bh_ffr_dm_run_query(BH_FFR_DM_PARTICLE_TYPE, dm);

  int local_need_fallback = 0;
  for(int n = 0; n < NumActiveBHFFR; n++)
    if(dm[n].Count < All.BHDMNeighbours)
      local_need_fallback = 1;

  int global_need_fallback = 0;
  MPI_Allreduce(&local_need_fallback, &global_need_fallback, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  struct bh_ffr_dm_result *gas = NULL;
  if(global_need_fallback)
    {
      gas = (struct bh_ffr_dm_result *)mymalloc(
          "BHFFRDMGasFallback", (NumActiveBHFFR > 0 ? NumActiveBHFFR : 1) * sizeof(*gas));
      bh_ffr_dm_run_query(BH_FFR_GAS_PARTICLE_TYPE, gas);
    }

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;

  for(int n = 0; n < NumActiveBHFFR; n++)
    {
      const int p = BHFFRActiveParticleList[n];
      const int b = P[p].BHDataIndex;

      const struct bh_ffr_dm_result *use = &dm[n];
      const char *source = "DM";

      if(dm[n].Count < All.BHDMNeighbours)
        {
          use = &gas[n];
          source = "gas-fallback";
        }

      const double sigma = bh_ffr_dm_sigma_1d(use);

      /* If even the fallback has fewer than two particles, retain an existing
       * valid cache (if any) rather than replacing it with a fabricated value.
       * A newly seeded BH therefore remains at SigmaDM=0 until a valid sample
       * is available; the optional central binding term still provides a
       * well-defined threshold in that situation. */
      if(use->Count >= 2)
        BHP[b].SigmaDM = sigma;

      /* Dynamical friction is specifically a DM wake model, so do not apply
       * it when SigmaDM had to fall back to gas. */
      if(use == &dm[n] && dm[n].Count >= 2)
        bh_ffr_apply_dynamical_friction(p, &dm[n]);

      if(!isfinite(BHP[b].SigmaDM) || BHP[b].SigmaDM < 0)
        terminate("BH_FFR: invalid cached SigmaDM=%g for ID=%llu", BHP[b].SigmaDM, (unsigned long long)P[p].ID);

      const double rmax_proper =
          use->Count > 0 ? sqrt(use->C[use->Count - 1].R2) * a : 0.0;

      printf("BH_FFR: sigmaDM cache ID=%llu task=%d source=%s count=%d requested=%d sigma=%g rmaxProper=%g\n",
             (unsigned long long)P[p].ID, ThisTask, source, use->Count, All.BHDMNeighbours, BHP[b].SigmaDM, rmax_proper);
      fflush(stdout);
    }

  if(gas != NULL)
    myfree(gas);
  myfree(dm);

  bh_ffr_validate_state("post-SigmaDM-cache");
}
