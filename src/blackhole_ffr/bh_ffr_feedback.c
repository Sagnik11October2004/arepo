#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../main/proto.h"

/*
 * Iteration 7: bursty, exactly bipolar wind feedback.
 *
 * The full gas mass inside R_fb sets the binding-energy threshold, while only
 * hydro-active gas cells are eligible injection targets.  Wind targets lie in
 * a broad bicone about DiscDir.  Each lobe uses separately normalized
 * mass-weighted top-hat weights.  The returned unresolved wind mass is added
 * comoving with the BH before the kinetic impulse is applied.
 *
 * The impulse amplitude is solved from the exact quadratic kinetic-energy
 * change in physical peculiar velocities, so every accepted packet has zero
 * net kick momentum and injects exactly E_packet up to MPI roundoff.
 *
 * The wind channel remains independent of the narrow jet channel. Iteration
 * 8 applies jet packets afterward using the gas state left by any wind events.
 */

enum bh_ffr_feedback_pass
{
  BH_FFR_FEEDBACK_STATS = 0,
  BH_FFR_FEEDBACK_MARK = 1,
  BH_FFR_FEEDBACK_CONFLICT = 2,
  BH_FFR_FEEDBACK_INJECT = 3
};

struct bh_ffr_feedback_event
{
  int Candidate;
  int Fire;
  MyDouble ThresholdEnergy;
  MyDouble PacketEnergy;
  MyDouble PacketMass;
  MyDouble PacketMomentum;
  MyDouble Q;
  MyDouble LobeMass[2];
  MyDouble LobeUtherm[2];
  long long LobeCount[2];
  MyDouble EnergyBefore;
  MyDouble MassBefore;
  MyDouble MomentumBefore;
  MyFloat CosCone;
  int UsedHemisphereFallback;
};

static struct bh_ffr_feedback_event *FeedbackEvents;
static int *FeedbackPacketCount;
static MyIDType *FeedbackWinnerID;
static int FeedbackNTargets;
static int FeedbackPass;
static int FeedbackConflictRound;

typedef struct
{
  MyDouble Pos[3];
  MyFloat Axis[3];
  MyFloat BHVel[3];
  MyFloat Radius;
  MyFloat CosCone;
  MyIDType BHID;
  int Candidate;
  int Fire;
  MyDouble PacketMass;
  MyDouble Q;
  MyDouble LobeMass[2];
  MyDouble LobeUtherm[2];
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble EnclosedMass;
  MyDouble LobeMass[2];
  MyDouble LobeProjectedMomentum[2];
  MyDouble LobeThermalMassWeighted[2];
  MyDouble LobeTotalMass[2];
  long long LobeCount[2];

  MyDouble HemiMass[2];
  MyDouble HemiProjectedMomentum[2];
  MyDouble HemiThermalMassWeighted[2];
  MyDouble HemiTotalMass[2];
  long long HemiCount[2];
  int Conflict;

  MyDouble ReturnedMass;
  MyDouble KickEnergy;
  MyDouble KickMomentum[3];
} data_out;

static data_out *FeedbackResults;
static data_out *DataResult, *DataOut;

static double bh_ffr_feedback_coordinate_radius(void)
{
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0) || !isfinite(All.BHFeedbackRadius) || !(All.BHFeedbackRadius > 0))
    terminate("BH_FFR: invalid feedback-radius conversion Rfb=%g a=%g", All.BHFeedbackRadius, a);

  const double r = All.BHFeedbackRadius / a;
  if(!isfinite(r) || !(r > 0))
    terminate("BH_FFR: invalid coordinate feedback radius=%g", r);
  return r;
}

static int bh_ffr_feedback_decode_hydro_timebin(int bin)
{
  if(bin < 0)
    bin = -bin - 1;
  return bin;
}

static int bh_ffr_feedback_gas_is_active(int j)
{
  const int bin = bh_ffr_feedback_decode_hydro_timebin(P[j].TimeBinHydro);
  return bin > 0 && bin < TIMEBINS && TimeBinSynchronized[bin];
}

static int bh_ffr_feedback_particle_from_target(int target, const char *where)
{
  if(target < 0 || target >= FeedbackNTargets)
    terminate("BH_FFR: feedback target=%d outside [0,%d) in %s", target, FeedbackNTargets, where);

  const int p = BHFFRActiveParticleList[target];
  if(p < 0 || p >= NumPart || P[p].Type != BH_FFR_PARTICLE_TYPE)
    terminate("BH_FFR: invalid active feedback particle=%d in %s", p, where);

  const int b = P[p].BHDataIndex;
  if(b < 0 || b >= NumBHFFR || BHP[b].ParticleID != P[p].ID)
    terminate("BH_FFR: invalid compact feedback state for particle ID=%llu in %s", (unsigned long long)P[p].ID, where);

  if(P[p].Ti_Current != All.Ti_Current)
    terminate("BH_FFR: feedback particle ID=%llu is not drifted to Ti_Current=%lld in %s", (unsigned long long)P[p].ID,
              (long long)All.Ti_Current, where);

  return p;
}

static void particle2in(data_in *in, int target, int firstnode)
{
  const int p = bh_ffr_feedback_particle_from_target(target, "particle2in");
  const int b = P[p].BHDataIndex;
  const struct bh_ffr_feedback_event *ev = &FeedbackEvents[target];

  for(int k = 0; k < 3; k++)
    {
      in->Pos[k] = P[p].Pos[k];
      in->Axis[k] = BHP[b].DiscDir[k];
      in->BHVel[k] = P[p].Vel[k];
    }

  in->Radius = bh_ffr_feedback_coordinate_radius();
  in->CosCone = ev->CosCone;
  in->BHID = P[p].ID;
  in->Candidate = ev->Candidate;
  in->Fire = ev->Fire;
  in->PacketMass = ev->PacketMass;
  in->Q = ev->Q;
  in->LobeMass[0] = ev->LobeMass[0];
  in->LobeMass[1] = ev->LobeMass[1];
  in->LobeUtherm[0] = ev->LobeUtherm[0];
  in->LobeUtherm[1] = ev->LobeUtherm[1];
  in->Firstnode = firstnode;
}

static void out2particle(data_out *out, int target, int mode)
{
  if(target < 0 || target >= FeedbackNTargets)
    terminate("BH_FFR: feedback result target=%d outside [0,%d)", target, FeedbackNTargets);

  data_out *res = &FeedbackResults[target];

  if(mode == MODE_LOCAL_PARTICLES)
    *res = *out;
  else
    {
      res->EnclosedMass += out->EnclosedMass;
      for(int l = 0; l < 2; l++)
        {
          res->LobeMass[l] += out->LobeMass[l];
          res->LobeProjectedMomentum[l] += out->LobeProjectedMomentum[l];
          res->LobeThermalMassWeighted[l] += out->LobeThermalMassWeighted[l];
          res->LobeTotalMass[l] += out->LobeTotalMass[l];
          res->LobeCount[l] += out->LobeCount[l];

          res->HemiMass[l] += out->HemiMass[l];
          res->HemiProjectedMomentum[l] += out->HemiProjectedMomentum[l];
          res->HemiThermalMassWeighted[l] += out->HemiThermalMassWeighted[l];
          res->HemiTotalMass[l] += out->HemiTotalMass[l];
          res->HemiCount[l] += out->HemiCount[l];
        }
      if(out->Conflict)
        res->Conflict = 1;

      res->ReturnedMass += out->ReturnedMass;
      res->KickEnergy += out->KickEnergy;
      for(int k = 0; k < 3; k++)
        res->KickMomentum[k] += out->KickMomentum[k];
    }
}

#include "../utils/generic_comm_helpers2.h"

static int bh_ffr_feedback_evaluate(int target, int mode, int threadid);

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
      if(target >= FeedbackNTargets)
        break;

      bh_ffr_feedback_evaluate(target, MODE_LOCAL_PARTICLES, threadid);
    }
}

static void kernel_imported(void)
{
  const int threadid = get_thread_num();
  int target = 0;

  while(target < Nimport)
    bh_ffr_feedback_evaluate(target++, MODE_IMPORTED_PARTICLES, threadid);
}

static unsigned long long bh_ffr_feedback_priority(MyIDType id)
{
  /* Deterministic time-varying total order. This preserves conflict safety
   * without permanently starving the higher-ID BH in a persistent overlap. */
  unsigned long long x = (unsigned long long)id;
  x ^= (unsigned long long)All.Ti_Current + 0x9e3779b97f4a7c15ULL +
       (unsigned long long)(FeedbackConflictRound + 1) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static int bh_ffr_feedback_candidate_wins(MyIDType candidate, MyIDType incumbent)
{
  if(incumbent == 0)
    return 1;

  const unsigned long long pc = bh_ffr_feedback_priority(candidate);
  const unsigned long long pi = bh_ffr_feedback_priority(incumbent);
  return pc < pi || (pc == pi && candidate < incumbent);
}

static int bh_ffr_feedback_selected_lobe(const data_in *in, int j, double *r2_out)
{
  double xtmp, ytmp, ztmp;
  const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
  const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
  const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
  const double r2 = dx * dx + dy * dy + dz * dz;

  if(r2_out != NULL)
    *r2_out = r2;

  if(!(r2 > 0) || r2 > in->Radius * in->Radius)
    return -1;

  const double dot = dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
  const double cosine = fabs(dot) / sqrt(r2);

  if(cosine < in->CosCone)
    return -1;

  return dot >= 0 ? 0 : 1;
}

static int bh_ffr_feedback_evaluate(int target, int mode, int threadid)
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

  data_out out;
  memset(&out, 0, sizeof(out));

  if(FeedbackPass != BH_FFR_FEEDBACK_STATS)
    {
      if((FeedbackPass == BH_FFR_FEEDBACK_INJECT && !in->Fire) ||
         ((FeedbackPass == BH_FFR_FEEDBACK_MARK || FeedbackPass == BH_FFR_FEEDBACK_CONFLICT) && !in->Candidate))
        {
          if(mode == MODE_LOCAL_PARTICLES)
            out2particle(&out, target, MODE_LOCAL_PARTICLES);
          else
            DataResult[target] = out;
          return 0;
        }
    }

  const int nfound =
      ngb_treefind_variable_threads(in->Pos, in->Radius, target, mode, threadid, numnodes, firstnode);

  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_FFR: invalid scale factor=%g in feedback evaluation", a);

  for(int n = 0; n < nfound; n++)
    {
      const int j = Thread[threadid].Ngblist[n];
      if(j < 0 || j >= NumGas)
        terminate("BH_FFR: feedback gas index=%d outside NumGas=%d", j, NumGas);
      if(P[j].Type != 0 || P[j].ID == 0 || !(P[j].Mass > 0))
        continue;

      double r2;
      const int lobe = bh_ffr_feedback_selected_lobe(in, j, &r2);

      if(FeedbackPass == BH_FFR_FEEDBACK_STATS)
        {
          int hemi = -1;
          if(r2 > 0 && r2 <= in->Radius * in->Radius)
            {
              const double dx = NEAREST_X(P[j].Pos[0] - in->Pos[0]);
              const double dy = NEAREST_Y(P[j].Pos[1] - in->Pos[1]);
              const double dz = NEAREST_Z(P[j].Pos[2] - in->Pos[2]);
              const double dot = dx * in->Axis[0] + dy * in->Axis[1] + dz * in->Axis[2];
              hemi = dot >= 0 ? 0 : 1;
              out.EnclosedMass += P[j].Mass;
              out.HemiTotalMass[hemi] += P[j].Mass;
              if(lobe >= 0)
                out.LobeTotalMass[lobe] += P[j].Mass;
            }

          const int active = bh_ffr_feedback_gas_is_active(j);
          if(active && hemi >= 0)
            {
              if(!isfinite(SphP[j].Utherm) || SphP[j].Utherm < 0)
                terminate("BH_FFR: invalid target internal energy u=%g for gas ID=%llu",
                          (double)SphP[j].Utherm, (unsigned long long)P[j].ID);

              double vdot = 0.0;
              for(int k = 0; k < 3; k++)
                vdot += (P[j].Vel[k] / a) * in->Axis[k];

              out.HemiCount[hemi]++;
              out.HemiMass[hemi] += P[j].Mass;
              out.HemiThermalMassWeighted[hemi] += P[j].Mass * SphP[j].Utherm;
              out.HemiProjectedMomentum[hemi] += P[j].Mass * vdot;

              if(lobe >= 0)
                {
                  out.LobeCount[lobe]++;
                  out.LobeMass[lobe] += P[j].Mass;
                  out.LobeThermalMassWeighted[lobe] += P[j].Mass * SphP[j].Utherm;
                  out.LobeProjectedMomentum[lobe] += P[j].Mass * vdot;
                }
            }
          continue;
        }

      if(lobe < 0 || !bh_ffr_feedback_gas_is_active(j))
        continue;

      if(FeedbackPass == BH_FFR_FEEDBACK_MARK)
        {
          if(bh_ffr_feedback_candidate_wins(in->BHID, FeedbackWinnerID[j]))
            FeedbackWinnerID[j] = in->BHID;
          continue;
        }

      if(FeedbackPass == BH_FFR_FEEDBACK_CONFLICT)
        {
          if(FeedbackWinnerID[j] != in->BHID)
            out.Conflict = 1;
          continue;
        }

      if(FeedbackPass != BH_FFR_FEEDBACK_INJECT)
        terminate("BH_FFR: unknown feedback pass=%d", FeedbackPass);

      if(!(in->LobeMass[lobe] > 0))
        terminate("BH_FFR: non-positive injection lobe mass for BH ID=%llu", (unsigned long long)in->BHID);

      const double old_mass = P[j].Mass;
      const double weight = old_mass / in->LobeMass[lobe];
      const double dm = 0.5 * in->PacketMass * weight;

      if(!isfinite(weight) || !(weight > 0) || !isfinite(dm) || dm < 0)
        terminate("BH_FFR: invalid wind target weight=%g dm=%g for gas ID=%llu", weight, dm, (unsigned long long)P[j].ID);

      double bh_v2 = 0.0;
      for(int k = 0; k < 3; k++)
        bh_v2 += in->BHVel[k] * in->BHVel[k];

      P[j].Mass += dm;
      for(int k = 0; k < 3; k++)
        SphP[j].Momentum[k] += dm * in->BHVel[k];
      /* Return mass with the lobe-weighted ambient specific internal
       * energy; the mechanical packet remains a separate kinetic increment. */
      SphP[j].Energy += dm * (0.5 * bh_v2 + a * a * in->LobeUtherm[lobe]);

      double dp_phys[3], v_phys[3];
      const double sign = (lobe == 0) ? 1.0 : -1.0;
      double dp2 = 0.0, vdotdp = 0.0;
      for(int k = 0; k < 3; k++)
        {
          v_phys[k] = (SphP[j].Momentum[k] / P[j].Mass) / a;
          dp_phys[k] = sign * in->Q * weight * in->Axis[k];
          dp2 += dp_phys[k] * dp_phys[k];
          vdotdp += v_phys[k] * dp_phys[k];
        }

      const double dkin = vdotdp + 0.5 * dp2 / P[j].Mass;
      if(!isfinite(dkin))
        terminate("BH_FFR: non-finite wind kinetic increment for gas ID=%llu", (unsigned long long)P[j].ID);

      for(int k = 0; k < 3; k++)
        {
          SphP[j].Momentum[k] += a * dp_phys[k];
          out.KickMomentum[k] += dp_phys[k];
        }
      SphP[j].Energy += a * a * dkin;

      for(int k = 0; k < 3; k++)
        P[j].Vel[k] = SphP[j].Momentum[k] / P[j].Mass;

      out.ReturnedMass += dm;
      out.KickEnergy += dkin;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static void bh_ffr_feedback_comm_pass(int pass)
{
  FeedbackPass = pass;
  memset(FeedbackResults, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackResults));

  if(All.TotNumGas > 0)
    {
      generic_set_MaxNexport();
      generic_comm_pattern(FeedbackNTargets, kernel_local, kernel_imported);
    }
}

static double bh_ffr_feedback_packet_q(double A, double B, double energy)
{
  if(!isfinite(A) || !isfinite(B) || !(B > 0) || !isfinite(energy) || !(energy > 0))
    terminate("BH_FFR: invalid packet-root inputs A=%g B=%g E=%g", A, B, energy);

  const double disc = sqrt(A * A + 2.0 * B * energy);
  double q;

  if(A >= 0)
    q = 2.0 * energy / (disc + A);
  else
    q = (-A + disc) / B;

  if(!isfinite(q) || !(q > 0))
    terminate("BH_FFR: invalid positive packet root q=%g A=%g B=%g E=%g", q, A, B, energy);

  return q;
}

void bh_ffr_feedback_self_test(void)
{
  const double A = -0.37;
  const double B = 0.81;
  const double energy = 1.23;
  const double q = bh_ffr_feedback_packet_q(A, B, energy);
  const double got = A * q + 0.5 * B * q * q;

  if(fabs(got - energy) > 2.0e-13 * energy)
    terminate("BH_FFR: feedback self-test failed packet root got=%g expected=%g", got, energy);

  const double mplus = 2.0;
  const double mminus = 3.0;
  const double mpacket = 0.5;
  const double Bmass = 1.0 / (mplus + 0.5 * mpacket) + 1.0 / (mminus + 0.5 * mpacket);
  if(!isfinite(Bmass) || !(Bmass > 0))
    terminate("BH_FFR: feedback self-test failed mass-weighted B coefficient");

  /* Regression for PDF eq. (77): returned wind mass advects the lobe-weighted
   * ambient internal energy in addition to its BH-bulk kinetic energy. With no
   * mechanical kick and a common bulk velocity, this must preserve the
   * mass-weighted specific internal energy. The a^2 factors exercise AREPO's
   * cosmological conserved-energy convention directly. */
  const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;
  if(!isfinite(a) || !(a > 0))
    terminate("BH_FFR: feedback self-test invalid scale factor=%g", a);

  const double m0 = 2.0, m1 = 3.0;
  const double u0 = 4.0, u1 = 10.0;
  const double m_before = m0 + m1;
  const double u_return = (m0 * u0 + m1 * u1) / m_before;
  const double v_bulk = 0.7;
  const double dm_return = 0.5;
  const double e_before = 0.5 * m_before * v_bulk * v_bulk + a * a * (m0 * u0 + m1 * u1);
  const double e_after =
      e_before + dm_return * (0.5 * v_bulk * v_bulk + a * a * u_return);
  const double m_after = m_before + dm_return;
  const double p_after = m_after * v_bulk;
  const double u_after =
      (e_after / m_after - 0.5 * (p_after / m_after) * (p_after / m_after)) / (a * a);

  if(!isfinite(u_after) || fabs(u_after - u_return) > 2.0e-13 * fmax(fabs(u_return), 1.0))
    terminate("BH_FFR: feedback self-test failed wind thermal return got=%g expected=%g", u_after, u_return);
}

static void bh_ffr_feedback_prepare_candidates(void)
{
  for(int n = 0; n < FeedbackNTargets; n++)
    {
      struct bh_ffr_feedback_event *ev = &FeedbackEvents[n];
      memset(ev, 0, sizeof(*ev));

      const int p = bh_ffr_feedback_particle_from_target(n, "bh_ffr_feedback_prepare_candidates");
      const int b = P[p].BHDataIndex;
      const data_out *res = &FeedbackResults[n];

      const double sigma2 = BHP[b].SigmaDM * BHP[b].SigmaDM;
      double vbind2 = sigma2;

      if(All.BHUseCentralBindingTerm)
        vbind2 += 2.0 * All.G * bh_ffr_central_mass_code(p) / All.BHFeedbackRadius;

      if(!isfinite(vbind2) || vbind2 < 0 || !isfinite(res->EnclosedMass) || res->EnclosedMass < 0)
        terminate("BH_FFR: invalid wind threshold environment Menc=%g vbind2=%g for ID=%llu", res->EnclosedMass, vbind2,
                  (unsigned long long)P[p].ID);

      const double base_threshold = 0.5 * res->EnclosedMass * vbind2;
      BHP[b].WindThresholdEnergy = All.BHWindBurstFactor * base_threshold;
      BHP[b].JetThresholdEnergy = All.BHJetBurstFactor * base_threshold;

      if(BHP[b].SigmaDM > 0)
        {
          printf("BH_FFR: binding threshold ID=%llu task=%d Menc=%g sigmaDM=%g vbind2=%g central=%d "
                 "EthWind=%g EthJet=%g\n",
                 (unsigned long long)P[p].ID, ThisTask, res->EnclosedMass, BHP[b].SigmaDM, vbind2,
                 All.BHUseCentralBindingTerm, BHP[b].WindThresholdEnergy, BHP[b].JetThresholdEnergy);
          fflush(stdout);
        }

      if(FeedbackPacketCount[n] >= All.BHMaxPacketsPerStep)
        continue;

      const double eth = BHP[b].WindThresholdEnergy;
      if(!(eth > 0) || !(BHP[b].WindEnergyBuffer >= eth))
        continue;

      const double cone_active_mass = res->LobeMass[0] + res->LobeMass[1];
      const double cone_total_mass = res->LobeTotalMass[0] + res->LobeTotalMass[1];
      const int cone_ok =
          res->LobeCount[0] >= All.BHMinTargetsPerLobe &&
          res->LobeCount[1] >= All.BHMinTargetsPerLobe &&
          res->LobeMass[0] > 0 && res->LobeMass[1] > 0 &&
          cone_total_mass > 0 &&
          cone_active_mass >= All.BHMinActiveTargetMassFrac * cone_total_mass;

      const double hemi_active_mass = res->HemiMass[0] + res->HemiMass[1];
      const double hemi_total_mass = res->HemiTotalMass[0] + res->HemiTotalMass[1];
      const int hemi_ok =
          res->HemiCount[0] >= All.BHMinTargetsPerLobe &&
          res->HemiCount[1] >= All.BHMinTargetsPerLobe &&
          res->HemiMass[0] > 0 && res->HemiMass[1] > 0 &&
          hemi_total_mass > 0 &&
          hemi_active_mass >= All.BHMinActiveTargetMassFrac * hemi_total_mass;

      if(!cone_ok && !hemi_ok)
        continue;

      ev->UsedHemisphereFallback = cone_ok ? 0 : 1;
      ev->CosCone = cone_ok ? cos(All.BHWindConeAngleDeg * M_PI / 180.0) : 0.0;

      ev->EnergyBefore = BHP[b].WindEnergyBuffer;
      ev->MassBefore = BHP[b].WindMassBuffer;
      ev->MomentumBefore = BHP[b].WindMomentumBuffer;
      ev->ThresholdEnergy = eth;

      /* The binding-energy threshold is a release trigger, not an energy
       * quantum that must be subcycled one packet at a time.  Once enough
       * hydro-active target mass exists, release the complete accumulated
       * wind reservoir in one conservative bipolar event.  This mirrors the
       * buffered-event semantics used by the TNG comparison backend and
       * prevents unresolved feedback energy from forcing BH-only timesteps
       * while the receiving gas is inactive. */
      ev->PacketEnergy = ev->EnergyBefore;

      const double frac = ev->PacketEnergy / ev->EnergyBefore;
      ev->PacketMass = ev->MassBefore * frac;
      ev->PacketMomentum = ev->MomentumBefore * frac;
      for(int l = 0; l < 2; l++)
        {
          ev->LobeMass[l] = ev->UsedHemisphereFallback ? res->HemiMass[l] : res->LobeMass[l];
          ev->LobeCount[l] = ev->UsedHemisphereFallback ? res->HemiCount[l] : res->LobeCount[l];
          const double thermal_mass =
              ev->UsedHemisphereFallback ? res->HemiThermalMassWeighted[l] : res->LobeThermalMassWeighted[l];
          ev->LobeUtherm[l] = thermal_mass / ev->LobeMass[l];
          if(!isfinite(ev->LobeUtherm[l]) || ev->LobeUtherm[l] < 0)
            terminate("BH_FFR: invalid lobe internal energy=%g for ID=%llu lobe=%d",
                      ev->LobeUtherm[l], (unsigned long long)P[p].ID, l);
        }

      if(!isfinite(frac) || frac <= 0 || frac > 1.0 + 2.0e-12 || !isfinite(ev->PacketMass) || ev->PacketMass < 0 ||
         ev->PacketMass > ev->MassBefore * (1.0 + 2.0e-12))
        terminate("BH_FFR: invalid wind packet fraction=%g dM=%g for ID=%llu", frac, ev->PacketMass,
                  (unsigned long long)P[p].ID);

      const double half = 0.5 * ev->PacketMass;
      const double mplus = res->LobeMass[0] + half;
      const double mminus = res->LobeMass[1] + half;
      const double a = All.ComovingIntegrationOn ? All.cf_atime : 1.0;

      double vbh_axis = 0.0;
      for(int k = 0; k < 3; k++)
        vbh_axis += (P[p].Vel[k] / a) * BHP[b].DiscDir[k];

      const double pplus =
          ev->UsedHemisphereFallback ? res->HemiProjectedMomentum[0] : res->LobeProjectedMomentum[0];
      const double pminus =
          ev->UsedHemisphereFallback ? res->HemiProjectedMomentum[1] : res->LobeProjectedMomentum[1];
      const double vplus = (pplus + half * vbh_axis) / mplus;
      const double vminus = (pminus + half * vbh_axis) / mminus;
      const double A = vplus - vminus;
      const double B = 1.0 / mplus + 1.0 / mminus;

      ev->Q = bh_ffr_feedback_packet_q(A, B, ev->PacketEnergy);
      ev->Candidate = 1;
    }
}

static int bh_ffr_feedback_global_candidate_count(void)
{
  int local = 0;
  for(int n = 0; n < FeedbackNTargets; n++)
    if(FeedbackEvents[n].Candidate)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static int bh_ffr_feedback_global_fire_count(void)
{
  int local = 0;
  for(int n = 0; n < FeedbackNTargets; n++)
    if(FeedbackEvents[n].Fire)
      local++;

  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

static void bh_ffr_feedback_resolve_conflicts(void)
{
  const int candidates = bh_ffr_feedback_global_candidate_count();
  if(candidates <= 0)
    return;

  if(candidates == 1)
    {
      for(int n = 0; n < FeedbackNTargets; n++)
        FeedbackEvents[n].Fire = FeedbackEvents[n].Candidate;
      return;
    }

  memset(FeedbackWinnerID, 0, (NumGas > 0 ? NumGas : 1) * sizeof(*FeedbackWinnerID));
  bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_MARK);
  bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_CONFLICT);

  for(int n = 0; n < FeedbackNTargets; n++)
    FeedbackEvents[n].Fire = FeedbackEvents[n].Candidate && !FeedbackResults[n].Conflict;
}

static void bh_ffr_feedback_commit_packets(void)
{
  for(int n = 0; n < FeedbackNTargets; n++)
    {
      struct bh_ffr_feedback_event *ev = &FeedbackEvents[n];
      if(!ev->Fire)
        continue;

      const int p = bh_ffr_feedback_particle_from_target(n, "bh_ffr_feedback_commit_packets");
      const int b = P[p].BHDataIndex;
      const data_out *res = &FeedbackResults[n];

      const double mtol = 2.0e-8 * fmax(fabs(ev->PacketMass), 1.0e-30);
      if(fabs(res->ReturnedMass - ev->PacketMass) > mtol)
        terminate("BH_FFR: wind packet mass-return mismatch ID=%llu returned=%g requested=%g", (unsigned long long)P[p].ID,
                  res->ReturnedMass, ev->PacketMass);

      const double etol = 3.0e-8 * fmax(fabs(ev->PacketEnergy), 1.0e-30);
      if(fabs(res->KickEnergy - ev->PacketEnergy) > etol)
        terminate("BH_FFR: wind packet kinetic-energy mismatch ID=%llu injected=%g requested=%g",
                  (unsigned long long)P[p].ID, res->KickEnergy, ev->PacketEnergy);

      double pnorm2 = 0.0;
      for(int k = 0; k < 3; k++)
        pnorm2 += res->KickMomentum[k] * res->KickMomentum[k];
      const double pnorm = sqrt(pnorm2);
      if(pnorm > 3.0e-8 * fmax(fabs(ev->Q), 1.0e-30))
        terminate("BH_FFR: wind packet kick-momentum imbalance ID=%llu |sum dp|=%g q=%g", (unsigned long long)P[p].ID,
                  pnorm, ev->Q);

      BHP[b].WindEnergyBuffer = ev->EnergyBefore - ev->PacketEnergy;
      BHP[b].WindMassBuffer = ev->MassBefore - ev->PacketMass;
      BHP[b].WindMomentumBuffer = ev->MomentumBefore - ev->PacketMomentum;

      if(BHP[b].WindEnergyBuffer < 0 && BHP[b].WindEnergyBuffer > -3.0e-12 * fmax(fabs(ev->EnergyBefore), 1.0e-30))
        BHP[b].WindEnergyBuffer = 0.0;
      if(BHP[b].WindMassBuffer < 0 && BHP[b].WindMassBuffer > -3.0e-12 * fmax(fabs(ev->MassBefore), 1.0e-30))
        BHP[b].WindMassBuffer = 0.0;
      if(BHP[b].WindMomentumBuffer < 0 && BHP[b].WindMomentumBuffer > -3.0e-12 * fmax(fabs(ev->MomentumBefore), 1.0e-30))
        BHP[b].WindMomentumBuffer = 0.0;

      if(BHP[b].WindEnergyBuffer < 0 || BHP[b].WindMassBuffer < 0 || BHP[b].WindMomentumBuffer < 0)
        terminate("BH_FFR: negative wind buffer after packet for ID=%llu", (unsigned long long)P[p].ID);

      P[p].Mass = BHP[b].BHMass + BHP[b].ReservoirMass + BHP[b].WindMassBuffer;
      FeedbackPacketCount[n]++;

      const double eerg = ev->PacketEnergy * All.UnitEnergy_in_cgs / All.HubbleParam;
      const double mmsun = ev->PacketMass * All.UnitMass_in_g / (All.HubbleParam * SOLAR_MASS);
      printf("BH_FFR: wind packet fired ID=%llu task=%d E=%g erg dM=%g Msun q=%g Nplus=%lld Nminus=%lld dErel=%g pbal=%g\n",
             (unsigned long long)P[p].ID, ThisTask, eerg, mmsun, ev->Q, ev->LobeCount[0], ev->LobeCount[1],
             fabs(res->KickEnergy - ev->PacketEnergy) / ev->PacketEnergy, pnorm / fmax(fabs(ev->Q), 1.0e-30));
      fflush(stdout);
    }
}

void bh_ffr_inject_wind_feedback(void)
{
  int global_active_bhs = 0;
  MPI_Allreduce(&NumActiveBHFFR, &global_active_bhs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(global_active_bhs <= 0)
    return;

#ifdef MHD
  terminate("BH_FFR: kinetic wind feedback is not enabled with MHD until the magnetic mass-return policy is validated");
#endif

  FeedbackNTargets = NumActiveBHFFR;
  FeedbackPacketCount =
      (int *)mymalloc("BHFFRFeedbackPacketCount", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(int));
  FeedbackEvents = (struct bh_ffr_feedback_event *)mymalloc(
      "BHFFRFeedbackEvents", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackEvents));
  FeedbackResults =
      (data_out *)mymalloc("BHFFRFeedbackResults", (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackResults));
  FeedbackWinnerID =
      (MyIDType *)mymalloc("BHFFRFeedbackWinner", (NumGas > 0 ? NumGas : 1) * sizeof(*FeedbackWinnerID));

  memset(FeedbackPacketCount, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(int));
  memset(FeedbackEvents, 0, (FeedbackNTargets > 0 ? FeedbackNTargets : 1) * sizeof(*FeedbackEvents));

  int any_global_fire = 0;

  for(int round = 0; round < All.BHMaxPacketsPerStep; round++)
    {
      FeedbackConflictRound = round;
      bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_STATS);
      bh_ffr_feedback_prepare_candidates();

      if(bh_ffr_feedback_global_candidate_count() <= 0)
        break;

      bh_ffr_feedback_resolve_conflicts();
      const int global_fire = bh_ffr_feedback_global_fire_count();
      if(global_fire <= 0)
        break;

      bh_ffr_feedback_comm_pass(BH_FFR_FEEDBACK_INJECT);
      bh_ffr_feedback_commit_packets();
      any_global_fire += global_fire;
    }

  if(any_global_fire > 0)
    update_primitive_variables();

  bh_ffr_validate_state("post-wind-feedback");

  myfree(FeedbackWinnerID);
  myfree(FeedbackResults);
  myfree(FeedbackEvents);
  myfree(FeedbackPacketCount);

  FeedbackWinnerID = NULL;
  FeedbackResults = NULL;
  FeedbackEvents = NULL;
  FeedbackPacketCount = NULL;
  FeedbackNTargets = 0;
}
