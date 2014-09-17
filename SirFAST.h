/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University, Bilkent University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the names of the University of Washington, Simon Fraser University,
 *   nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/*
  Authors:
  Farhad Hormozdiari
  Faraz Hach
  Can Alkan
  Emails:
  farhadh AT uw DOT edu
  fhach AT cs DOT sfu DOT ca
  calkan AT cs DOT bilkent DOT edu DOT tr
*/



#ifndef __MR_FAST__
#define __MR_FAST__

#include "Reads.h"
#include <ctype.h>

#define MAP_CHUNKS 15
#define MAX_CIGAR_SIZE 100
#define MAX_CG_SIZE 20


// Pair is used to pre-processing and making the read index table
typedef struct
{
  int hv;
  //char hv[50];
  int readNumber;
} Pair;

typedef struct
{
  int hv;
  unsigned int *seqInfo;
} ReadIndexTable;


typedef struct
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
} FullMappingInfo;

typedef struct
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char chr[MAX_CIGAR_SIZE];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
  double tprob;
} BestFullMappingInfo;

typedef struct lc
{
  char md[MAP_CHUNKS][MAX_CG_SIZE];
  int mdSize[MAP_CHUNKS];

  char cigar[MAP_CHUNKS][MAX_CG_SIZE];
  int cigarSize[MAP_CHUNKS];

  int err[MAP_CHUNKS];
  int loc[MAP_CHUNKS];
  struct lc *next;
} MappingLocations;

typedef struct inf
{
  int size;
  MappingLocations *next;
} MappingInfo;


typedef struct
{
  FullMappingInfo *mi;
  int size;
} FullMappingInfoLink;


extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;

void initFAST(Read *, int, int *, int, char *);

void initVerifiedLocs();
void initLookUpTable();


void finalizeFAST();
void finalizeBestSingleMapping();
void finalizeBestConcordantDiscordant();
//void finalizeOEAReads(char *);


int mapAllSingleEndSeq();

void generateCigarFromMD(char *, int, char *);

int msfHashVal(char *);

int backwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);
int forwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);

double mapProb(int, char *, int, int);
int mapQ(int);

/***********************************/

// for fastHASH
int compareEntrySize (const void *a, const void *b);											// fastHASH()
void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,	// fastHASH()
					 int index, key_struct* keys_input, int potential_key_number); 				// fastHASH()
void mapPairEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,// fastHASH()
					   int index, key_struct* keys_input, int potential_key_number); 			// fastHASH(
void mapPairedEndSeq();
//void outputPairedEnd();
int outputPairedEnd(int pre_unmappedCnt);
void setFullMappingInfo(int readNumber, int loc, int dir, int err, int score,
			char *md, char * refName, char *cigar);

void outputAllTransChromosomal(int flag);
//void outputPairedEndDiscPP();

// sirFAST
void initFASTCG(Read *seqList, int seqListSize, int accSeqListSize);
void preProcessReadsCG();
void mapAllSingleEndSeqCG();
int searchKeyCG(int target_coor, unsigned int* entry_coor, int entry_size, int range);

void mapSingleEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext);
void mapSingleEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext);
void mapPairEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext);
void mapPairEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext);
void printoutCG_forward(int genLoc, int * af_offset, char * seq, int error);
void printoutCG_reverse(int genLoc, int * af_offset, char * seq, int error);

int verifySingleEndCG(int refIndex, int readIndex, char* seq, int * af_offset, int variable, int length);

int verifySingleEndSeqCG_backward(int * locs, int * af_offset, char * seq1, key_struct* key_input, int index, int num_base, int num_key, int num_ver, int num_ext);
int verifySingleEndSeqCG_forward(int * locs, int * af_offset, char * seq1, key_struct* key_input, int index, int num_base, int num_key, int num_ver, int num_ext);

void generateAlignmentMatrxCG_backward(int genLoc, int * af_offset, char * seq, int error, char * matrix);
void generateAlignmentMatrxCG_forward(int genLoc, int * af_offset, char * seq, int error, char * matrix);

void resetFAST(unsigned int seqListSize);

FILE * initOEAReads();
void finalizeOEAReads(FILE * _fp_oea);
void performOEAReads(FILE * _fp_oea);


FILE * initPairedEndDiscPP();
void finalizePairedEndDiscPP(FILE * _fp_divet);
void performPairedEndDiscPP(FILE * _fp_divet);

void initBestMapping();
void performBestMapping(int readCnt);
void finalizeBestMapping(int readCnt);
void mapPairedEndSeqCG();

void initKeyInput(key_struct* key_input, int* key_hash);
void initForwardParam();
void initReverseParam();

/**************************************************************************************************/

#endif
