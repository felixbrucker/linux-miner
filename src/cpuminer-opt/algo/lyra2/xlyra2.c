
int LYRA2REV2( uint64_t* wholeMatrix, void *K, uint64_t kLen, const void *pwd,
               uint64_t pwdlen, const void *salt, uint64_t saltlen,
               uint64_t timeCost, const uint64_t nRows, const uint64_t nCols )
{
   //====================== Basic variables ============================//
   uint64_t _ALIGN(256) state[16];
   int64_t row = 2; //index of row to be processed
   int64_t prev = 1; //index of prev (last row ever computed/modified)
   int64_t rowa = 0; //index of row* (a previous row, deterministically picked during Setup and randomly picked while Wandering)
   int64_t tau; //Time Loop iterator
   int64_t step = 1; //Visitation step (used during Setup and Wandering phases)
   int64_t window = 2; //Visitation window (used to define which rows can be revisited during Setup)
   int64_t gap = 1; //Modifier to the step, assuming the values 1 or -1
   int64_t i; //auxiliary iteration counter
   int64_t v64; // 64bit var for memcpy
   //====================================================================/

   //=== Initializing the Memory Matrix and pointers to it =============//
   //Tries to allocate enough space for the whole memory matrix

   const int64_t ROW_LEN_INT64 = BLOCK_LEN_INT64 * nCols;
   const int64_t ROW_LEN_BYTES = ROW_LEN_INT64 * 8;
   // for Lyra2REv2, nCols = 4, v1 was using 8
   const int64_t BLOCK_LEN = (nCols == 4) ? BLOCK_LEN_BLAKE2_SAFE_INT64
                                          : BLOCK_LEN_BLAKE2_SAFE_BYTES;
/*
   i = (int64_t)ROW_LEN_BYTES * nRows;
   uint64_t *wholeMatrix = _mm_malloc( i, 64 );
   if (wholeMatrix == NULL)
      return -1;

#if defined (__AVX2__)
   memset_zero_m256i( (__m256i*)wholeMatrix, i/32 );
#elif defined(__AVX__)
   memset_zero_m128i( (__m128i*)wholeMatrix, i/16 );
#else
   memset(wholeMatrix, 0, i);
#endif
*/
   uint64_t *ptrWord = wholeMatrix;

   //=== Getting the password + salt + basil padded with 10*1 ==========//
   //OBS.:The memory matrix will temporarily hold the password: not for saving memory,
   //but this ensures that the password copied locally will be overwritten as soon as possible

   //First, we clean enough blocks for the password, salt, basil and padding
   int64_t nBlocksInput = ( ( saltlen + pwdlen + 6 * sizeof(uint64_t) )
                              / BLOCK_LEN_BLAKE2_SAFE_BYTES ) + 1;

   byte *ptrByte = (byte*) wholeMatrix;

   //Prepends the password
   memcpy(ptrByte, pwd, pwdlen);
   ptrByte += pwdlen;

   //Concatenates the salt
   memcpy(ptrByte, salt, saltlen);
   ptrByte += saltlen;

   memset( ptrByte, 0, nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES
                       - (saltlen + pwdlen) );

   //Concatenates the basil: every integer passed as parameter, in the order they are provided by the interface
   memcpy(ptrByte, &kLen, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);
   v64 = pwdlen;
   memcpy(ptrByte, &v64, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);
   v64 = saltlen;
   memcpy(ptrByte, &v64, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);
   v64 = timeCost;
   memcpy(ptrByte, &v64, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);
   v64 = nRows;
   memcpy(ptrByte, &v64, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);
   v64 = nCols;
   memcpy(ptrByte, &v64, sizeof(int64_t));
   ptrByte += sizeof(uint64_t);

   //Now comes the padding
   *ptrByte = 0x80; //first byte of padding: right after the password
   ptrByte = (byte*) wholeMatrix; //resets the pointer to the start of the memory matrix
   ptrByte += nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - 1; //sets the pointer to the correct position: end of incomplete block
   *ptrByte ^= 0x01; //last byte of padding: at the end of the last incomplete block

// from here on it's all simd acces to state and matrix
// define vector pointers and adjust sizes and pointer offsets

   //================= Initializing the Sponge State ====================//
   //Sponge state: 16 uint64_t, BLOCK_LEN_INT64 words of them for the bitrate (b) and the remainder for the capacity (c)

   initState( state );

   //========================= Setup Phase =============================//
   //Absorbing salt, password and basil: this is the only place in which the block length is hard-coded to 512 bits
   
   ptrWord = wholeMatrix;
   for (i = 0; i < nBlocksInput; i++)
   {
       absorbBlockBlake2Safe( state, ptrWord ); //absorbs each block of pad(pwd || salt || basil)
       ptrWord += BLOCK_LEN; //goes to next block of pad(pwd || salt || basil)
   }
   //Initializes M[0] and M[1]
   reducedSqueezeRow0( state, &wholeMatrix[0], nCols ); //The locally copied password is most likely overwritten here

   reducedDuplexRow1( state, &wholeMatrix[0], &wholeMatrix[ROW_LEN_INT64],
                      nCols);

   do
   {
      //M[row] = rand; //M[row*] = M[row*] XOR rotW(rand)

      reducedDuplexRowSetup( state, &wholeMatrix[prev*ROW_LEN_INT64],
                             &wholeMatrix[rowa*ROW_LEN_INT64],
                             &wholeMatrix[row*ROW_LEN_INT64], nCols );

      //updates the value of row* (deterministically picked during Setup))
      rowa = (rowa + step) & (window - 1);
      //update prev: it now points to the last row ever computed

      prev = row;
      //updates row: goes to the next row to be computed
      row++;

      //Checks if all rows in the window where visited.
      if (rowa == 0)
      {
         step = window + gap; //changes the step: approximately doubles its value
         window *= 2; //doubles the size of the re-visitation window
         gap = -gap; //inverts the modifier to the step
      }

   } while (row < nRows);

   //===================== Wandering Phase =============================//
   row = 0; //Resets the visitation to the first row of the memory matrix
   for (tau = 1; tau <= timeCost; tau++)
   {
       //Step is approximately half the number of all rows of the memory matrix for an odd tau; otherwise, it is -1
       step = (tau % 2 == 0) ? -1 : nRows / 2 - 1;
       do
       {
           //Selects a pseudorandom index row*
           //-----------------------------------------------
           rowa = state[0] & (unsigned int)(nRows-1);  //(USE THIS IF nRows IS A POWER OF 2)

           //rowa = state[0] % nRows; //(USE THIS FOR THE "GENERIC" CASE)
           //-------------------------------------------

           //Performs a reduced-round duplexing operation over M[row*] XOR M[prev], updating both M[row*] and M[row]
           reducedDuplexRow( state, &wholeMatrix[prev*ROW_LEN_INT64],
                             &wholeMatrix[rowa*ROW_LEN_INT64],
                             &wholeMatrix[row*ROW_LEN_INT64], nCols );
           //update prev: it now points to the last row ever computed
           prev = row;

           //updates row: goes to the next row to be computed
           //----------------------------------------------------
           row = (row + step) & (unsigned int)(nRows-1); //(USE THIS IF nRows IS A POWER OF 2)
           //row = (row + step) % nRows; //(USE THIS FOR THE "GENERIC" CASE)
           //----------------------------------------------------

       } while (row != 0);
   }

   //===================== Wrap-up Phase ===============================//
   //Absorbs the last block of the memory matrix
   absorbBlock(state, &wholeMatrix[rowa*ROW_LEN_INT64]);
   //Squeezes the key
   squeeze(state, K, (unsigned int) kLen);

   //================== Freeing the memory =============================//
//   free(wholeMatrix);

   return 0;
}

