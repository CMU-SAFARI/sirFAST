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

int compare(const void *a, const void *b); 
float str2int(char *str, int index1, int index2);
void initBestMapping(int totalReadNumber);
void finalizeBestSingleMapping() ;
void preProcessReads();
void resetSingleFAST(unsigned int seqListSize);
void resetFAST(unsigned int seqListSize);
void finalizeFAST();
int addCigarSize(int cnt);
void generateSNPSAM(char *matrix, int matrixLength, char *outputSNP);
int compareOut(const void *a, const void *b);
void outputPairFullMappingInfo(int readNumber);
int findNearest(int x1, int x2, int c);
void finalizeBestConcordantDiscordant();
double mapProb(int readNumber, char *md, int dir, int err);
int mapQ(int readNumber);
void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1, FullMappingInfo mi2);
int outputPairedEnd(int pre_unmappedCnt);
float calculateScore(int index, char *seq, char *qual, char *md);
int matoi(char *str, int start, int end); 
void convertCigarToMatrix(char *cigar, int cigar_size, char * matrix);
void convertMDToMatrix(char *md, int md_size, char * matrix);
void convertMDCigarToMatrix(char *cigar, int cigar_size, char *md, int md_size, char *matrix);
void convertInsertion(char * in_matrix, char * seq, char *out_matrix);
void finalizePairedEndDiscPP(FILE * out); 
void operatePairedEndDiscPP(FILE * out); 
void finalizeOEAReads(FILE * fp_out1); 
void operateOEAReads(FILE * fp_out1); 
void outputAllTransChromosomal(int flag);
void initFASTCG(Read *seqList, int seqListSize, char *genFileName, int AccReads, int first_try);
void mapPairedEndSeqCG(Read *seqList, unsigned int seqListSize, unsigned int AccReads); 
int searchKeyCG(int target_coor, unsigned int* entry_coor, int entry_size, int range); 
int verifySingleEndCG1(int refIndex, char* seq1, int * tmp_offset, int variable, int length);
int verifySingleEndCG(int refIndex, int readIndex, char * seq1, int * offset, int variable, int length);
void mapAllSingleEndSeqCG();
void freeAllMapping();
void mapPairEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int first_mate, int thread_id);
void mapPairEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int first_mate, int thread_id);
void generateAlignmentMatrxCG_backward(int genLoc, int * af_offset, char * seq, int * af_pass, int error, char * matrix);
void generateAlignmentMatrxCG_forward(int genLoc, int * af_offset, char * seq, int * af_pass, int error, char * matrix);
int verifySingleEndSeqCG_forward(int * locs, int * af_offset, char * seq1, int * af_pass, key_struct* key_input, int index);
int verifySingleEndSeqCG_backward(int * locs, int * af_offset, char * seq1, int * af_pass, key_struct* key_input, int index);
void mapSingleEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int first_mate, int thread_id);
void mapSingleEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index, key_struct* key_input, int direction, int first_mate, int thread_id);
FILE * initOEAReads(char *fileName);
FILE * initPairedEndDiscPP();

#endif
