#include <stdlib.h>
#include <string.h>

#include "blackhole_ffr.h"
#include "../domain/domain.h"
#include "../main/proto.h"

static int bh_ffr_domain_target_for_particle(int i)
{
  int no = 0;
  peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
  int shift = 3 * (BITS_PER_DIMENSION - 1);

  while(topNodes[no].Daughter >= 0)
    {
      no = topNodes[no].Daughter + (int)((Key[i] & mask) >> shift);
      mask >>= 3;
      shift -= 3;
    }

  return DomainTask[topNodes[no].Leaf];
}

void bh_ffr_domain_exchange_begin(struct bh_ffr_domain_exchange_context *ctx)
{
  if(ctx == NULL)
    terminate("BH_FFR: NULL domain-exchange context");

  ctx->Received = NULL;
  ctx->NumReceived = 0;

  bh_ffr_validate_state("domain-exchange begin");

  int *send_count = (int *)mymalloc("BHFFRSendCount", NTask * sizeof(int));
  int *recv_count = (int *)mymalloc("BHFFRRecvCount", NTask * sizeof(int));
  int *filled = (int *)mymalloc("BHFFRFilled", NTask * sizeof(int));
  size_t *send_offset = (size_t *)mymalloc("BHFFRSendOffset", NTask * sizeof(size_t));
  size_t *recv_offset = (size_t *)mymalloc("BHFFRRecvOffset", NTask * sizeof(size_t));
  size_t *send_count_size = (size_t *)mymalloc("BHFFRSendCountSize", NTask * sizeof(size_t));
  size_t *recv_count_size = (size_t *)mymalloc("BHFFRRecvCountSize", NTask * sizeof(size_t));

  memset(send_count, 0, NTask * sizeof(int));
  memset(recv_count, 0, NTask * sizeof(int));
  memset(filled, 0, NTask * sizeof(int));

  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE)
      {
        const int target = bh_ffr_domain_target_for_particle(i);
        if(target != ThisTask)
          send_count[target]++;
      }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  size_t nsend = 0;
  size_t nrecv = 0;
  for(int task = 0; task < NTask; task++)
    {
      send_offset[task] = nsend;
      recv_offset[task] = nrecv;
      send_count_size[task] = (size_t)send_count[task];
      recv_count_size[task] = (size_t)recv_count[task];
      nsend += send_count_size[task];
      nrecv += recv_count_size[task];
    }

  struct bh_ffr_particle_data *sendbuf =
      (struct bh_ffr_particle_data *)mymalloc("BHFFRSendBuf", (nsend > 0 ? nsend : 1) * sizeof(struct bh_ffr_particle_data));
  /* This receive buffer has to remain alive while AREPO's domain exchange
   * allocates and frees its own movable buffers. It therefore cannot live in
   * the non-movable mymalloc stack. */
  ctx->Received =
      (struct bh_ffr_particle_data *)malloc((nrecv > 0 ? nrecv : 1) * sizeof(struct bh_ffr_particle_data));
  if(ctx->Received == NULL)
    terminate("BH_FFR: failed to allocate domain receive buffer for %zu records", nrecv);
  ctx->NumReceived = (int)nrecv;

  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE)
      {
        const int target = bh_ffr_domain_target_for_particle(i);
        if(target != ThisTask)
          {
            const int b = P[i].BHDataIndex;
            const size_t slot = send_offset[target] + (size_t)filled[target]++;
            sendbuf[slot] = BHP[b];
          }
      }

  myMPI_Alltoallv(sendbuf, send_count_size, send_offset, ctx->Received, recv_count_size, recv_offset,
                  sizeof(struct bh_ffr_particle_data), 0, MPI_COMM_WORLD);

  myfree(sendbuf);
  myfree(recv_count_size);
  myfree(send_count_size);
  myfree(recv_offset);
  myfree(send_offset);
  myfree(filled);
  myfree(recv_count);
  myfree(send_count);
}

void bh_ffr_domain_exchange_finish(struct bh_ffr_domain_exchange_context *ctx)
{
  if(ctx == NULL || ctx->Received == NULL)
    terminate("BH_FFR: invalid domain-exchange context at finish");

  struct bh_ffr_particle_data *old_data = BHP;
  const int old_count = NumBHFFR;

  int new_count = 0;
  for(int i = 0; i < NumPart; i++)
    if(P[i].Type == BH_FFR_PARTICLE_TYPE)
      new_count++;

  struct bh_ffr_particle_data *new_data = NULL;
  if(new_count > 0)
    {
      new_data = (struct bh_ffr_particle_data *)malloc(new_count * sizeof(struct bh_ffr_particle_data));
      if(new_data == NULL)
        terminate("BH_FFR: failed to allocate post-domain compact state for %d black holes", new_count);
    }

  unsigned char *old_used = NULL;
  unsigned char *recv_used = NULL;

  if(old_count > 0)
    {
      old_used = (unsigned char *)mymalloc("BHFFROldUsed", old_count * sizeof(unsigned char));
      memset(old_used, 0, old_count * sizeof(unsigned char));
    }

  if(ctx->NumReceived > 0)
    {
      recv_used = (unsigned char *)mymalloc("BHFFRRecvUsed", ctx->NumReceived * sizeof(unsigned char));
      memset(recv_used, 0, ctx->NumReceived * sizeof(unsigned char));
    }

  int bnew = 0;
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].Type != BH_FFR_PARTICLE_TYPE)
        {
          P[i].BHDataIndex = -1;
          continue;
        }

      int found = 0;
      const int old_index = P[i].BHDataIndex;

      if(old_index >= 0 && old_index < old_count && old_data[old_index].ParticleID == P[i].ID && !old_used[old_index])
        {
          new_data[bnew] = old_data[old_index];
          old_used[old_index] = 1;
          found = 1;
        }

      if(!found)
        for(int j = 0; j < ctx->NumReceived; j++)
          if(!recv_used[j] && ctx->Received[j].ParticleID == P[i].ID)
            {
              new_data[bnew] = ctx->Received[j];
              recv_used[j] = 1;
              found = 1;
              break;
            }

      if(!found)
        terminate("BH_FFR: no compact state found after domain exchange for particle ID=%llu", (unsigned long long)P[i].ID);

      P[i].BHDataIndex = bnew;
      bnew++;
    }

  for(int j = 0; j < ctx->NumReceived; j++)
    if(!recv_used[j])
      terminate("BH_FFR: received compact state for ID=%llu was not matched to a local Type-5 particle",
                (unsigned long long)ctx->Received[j].ParticleID);

  /* old_used and recv_used are AREPO-stack scratch blocks, so unwind them
   * in exact reverse allocation order. */
  if(recv_used != NULL)
    myfree(recv_used);
  if(old_used != NULL)
    myfree(old_used);

  if(old_data != NULL)
    free(old_data);

  BHP = new_data;
  NumBHFFR = new_count;

  free(ctx->Received);
  ctx->Received = NULL;
  ctx->NumReceived = 0;

  bh_ffr_validate_state("domain-exchange finish");
}
