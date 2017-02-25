/**
 * Blake2-B CUDA Implementation
 *
 * tpruvot@github July 2016
 *
 */

#include <miner.h>

#include <string.h>
#include <stdint.h>

#include <sph_blake2b.h>


#define NBN 2


void blake2b_hash(void *output, const void *input)
{
	uint8_t _ALIGN(A) hash[32];
	blake2b_ctx ctx;

	blake2b_init(&ctx, 32, NULL, 0);
	blake2b_update(&ctx, input, 80);
	blake2b_final(&ctx, hash);

	memcpy(output, hash, 32);
}

// ----------------------------------------------------------------


int scanhash_blake2b( int thr_id, struct work *work, uint32_t max_nonce,
                  uint32_t *hashes_done )
{
	uint32_t _ALIGN(64) hash[8];
	uint32_t _ALIGN(64) inputdata[20];
	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;

	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[8];
        uint32_t n = first_nonce;

	memcpy(inputdata, pdata, 80);
	inputdata[11] = 0; // nbits

	do {
		blake2b_hash(hash, inputdata);
		if (swab32(hash[0]) <= Htarg)
                {
			// sia hash target is reversed (start of hash)
//			swab256(vhashcpu, hash);
			if (fulltest(hash, ptarget))
                        {
                                *hashes_done = n - first_nonce + 1;
                                return true;
			}
		}
                n++; pdata[8] = n;

	} while (!work_restart[thr_id].restart);

	*hashes_done = pdata[8] - first_nonce;

	return 0;
}

/* compute nbits to get the network diff */
void blake2b_calc_network_diff(struct work *work)
{
        // sample for diff 43.281 : 1c05ea29
        // todo: endian reversed on longpoll could be zr5 specific...
        uint32_t nbits = have_longpoll ? work->data[18] : swab32(work->data[18]);
        if (opt_algo == ALGO_LBRY) nbits = swab32(work->data[26]);
        if (opt_algo == ALGO_DECRED) nbits = work->data[29];
        if (opt_algo == ALGO_SIA) nbits = work->data[11]; // unsure if correct

        uint32_t bits = (nbits & 0xffffff);
        int16_t shift = (swab32(nbits) & 0xff); // 0x1c = 28

        uint64_t diffone = 0x0000FFFF00000000ull;
        double d = (double)0x0000ffff / (double)bits;

        for (int m=shift; m < 29; m++) d *= 256.0;
        for (int m=29; m < shift; m++) d /= 256.0;
        if (opt_algo == ALGO_DECRED && shift == 28) d *= 256.0;
        if (opt_debug_diff)
                applog(LOG_DEBUG, "net diff: %f -> shift %u, bits %08x", d, shift, bits);

        net_diff = d;
}

void blake2b_le_build_stratum_request( char *req, struct work *work )
{
   unsigned char *xnonce2str;
   uint32_t ntime,       nonce;
   char     ntimestr[9], noncestr[9];
   le32enc( &ntime, work->data[ algo_gate.ntime_index ] );
   le32enc( &nonce, work->data[ algo_gate.nonce_index ] );
   bin2hex( ntimestr, (char*)(&ntime), sizeof(uint32_t) );
   bin2hex( noncestr, (char*)(&nonce), sizeof(uint32_t) );
   uint16_t high_nonce = swab32(work->data[9]) >> 16;
   xnonce2str = bin2hex((unsigned char*)(&high_nonce), 2);
//   xnonce2str = abin2hex(work->xnonce2, work->xnonce2_len);
   snprintf( req, JSON_BUF_LEN,
        "{\"method\": \"mining.submit\", \"params\": [\"%s\", \"%s\", \"%s\", \"%s\", \"%s\"], \"id\":4}",
         rpc_user, work->job_id, xnonce2str, ntimestr, noncestr );
   free( xnonce2str );
}

void blake2b_build_extraheader( struct work* work, struct stratum_ctx* sctx )
{
    uint32_t extra = 0;
    memcpy(&extra, &sctx->job.coinbase[32], 2);
    for (i = 0; i < 8; i++) // reversed hash
        work->data[i] = ((uint32_t*)sctx->job.prevhash)[7-i];
    work->data[8] = 0; // nonce
    work->data[9] = swab32(extra) | ((rand() << 8) & 0xffff);
    work->data[10] = be32dec(sctx->job.ntime);
    work->data[11] = be32dec(sctx->job.nbits);
    memcpy(&work->data[12], sctx->job.coinbase, 32); // merkle_root
    work->data[20] = 0x80000000;
}

void blake2b_get_new_work( struct work* work, struct work* g_work, int thr_id,
                           uint32_t* end_nonce_ptr, bool clean_job )
{
   const int wkcmp_sz = 32;
   const int wkcmp_off = ( 32 + 16 ) / 4;
   uint32_t *nonceptr = algo_gate.get_nonceptr( work->data );
   if ( memcmp( &work->data[ wkcmp_off ], &g_work->data[1], wkcmp_sz )
      && ( clean_job || ( *nonceptr >= *end_nonce_ptr ) ) )
      work_free( work );
      work_copy( work, g_work );
      *nonceptr = ( 0xffffffffU / opt_n_threads ) * thr_id;
      if ( opt_randomize )
         *nonceptr += ( (rand() *4 ) & UINT32_MAX ) / opt_n_threads;
      *end_nonce_ptr = ( 0xffffffffU / opt_n_threads ) * (thr_id+1) - 0x20;
   }
   else
       ++(*nonceptr);

   // suprnova job_id check without data/target/height change...
   if ( have_stratum && strcmp( work.job_id, g_work.job_id ) )
   {
       pthread_mutex_unlock(&g_work_lock);
       work_done = true;
       continue;
   }
   nonceptr[1] += opt_n_threads;
   nonceptr[1] |= thr_id;
   // range max
   nonceptr[0] = 0;
   end_nonce = UINT32_MAX;


bool register_blake2b_algo( algo_gate_t* gate )
{
  gate->ntime_index = 10;
  gate->nbits_index = 11;
  gate->nonce_index = 8;
  gate->scanhash              = (void*)&scanhash_blake2b;
  gate->hash                  = (void*)&blake2bhash;
  gate->calc_network_diff     = (void*)&blake2b_cal_network_diff;
  gate->build_stratum_request = (void*)&blake2b_le_build_stratum_request;
  gate->gen_merkle_root       = (void*)&do_nothing;
  gate->build_extraheader     = (void*)&blake2b_build_extraheader;
//  gate->hash_alt  = (void*)&blakehash;
//  gate->get_max64 = (void*)&blake_get_max64;
  return true;
}
