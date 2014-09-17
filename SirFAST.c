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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "Common.h"
#include "Reads.h"
#include "HashTable.h"
#include "Output.h"
#include "SirFAST.h"
#include "RefGenome.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define min3(a,b,c) ((a)>(b)?(b>c?c:b):(a>c?c:a))
#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))

#define MAX_REF_SIZE	 18
#define KEY_LENGTH     10              // sirfast
#define KEY_LENGTH0     5              // sirfast
#define INDEL_GAP       3              // sirfast
#define EXPECTED_GAP    2              // sirfast
#define INITIAL_GAP01  -1              // sirfast
#define INITIAL_GAP12   0              // sirfast
#define INITIAL_GAP23   5              // sirfast
#define DEBUG           3
#define TEST_KEY_NUM	 	3
#define READ_NUM	 			4

int key_length[4];
int initial_gap[3];
int expected_gap[3];

double binomial_coefficient(int n, int k);
float calculateScore(int index, char *seq, char *qual, char *md);

char *versionNumberF = "0.4";

long long verificationCnt = 0;
long long mappingCnt = 0;
long long mappedSeqCnt = 0;
long long completedSeqCnt = 0;
char *mappingOutput;
char *_msf_refGen = NULL;
int _msf_refGenLength = 0;
int _msf_refGenOffset = 0;
char *_msf_refGenName = NULL;
int _msf_refGenBeg;
int _msf_refGenEnd;
IHashTable *_msf_hashTable = NULL;
int *_msf_samplingLocs;
int *_msf_samplingLocsEnds;
int _msf_samplingLocsSize;
Read *_msf_seqList;
int _msf_seqListSize;
Pair *_msf_sort_seqList = NULL;
int *_msf_map_sort_seqList;
ReadIndexTable *_msf_rIndex = NULL;
int _msf_rIndexSize;
int _msf_rIndexMax;
SAM _msf_output;
OPT_FIELDS *_msf_optionalFields;
char *_msf_op;
int *_msf_verifiedLocs = NULL;
char _msf_numbers[200][3];
char _msf_cigar[5];
// READ SLICE
int _msf_totalSeqListSize;


MappingInfo *_msf_mappingInfo;

int *_msf_seqHits;
int _msf_openFiles = 0;
int _msf_maxLSize = 0;
int _msf_maxRSize = 0;

BestFullMappingInfo *bestHitMappingInfo;

int _msf_maxFile = 0;
char _msf_fileName[4000][200][2][FILE_NAME_LENGTH];
int _msf_fileCount[4000];

char *_msf_readHasConcordantMapping;

int *_msf_oeaMapping;
int *_msf_discordantMapping;

void freeAllMapping();

void initBestMapping(int totalReadNumber) {
	int i = 0;
	bestHitMappingInfo = getMem(totalReadNumber * sizeof(BestFullMappingInfo));
	for (i = 0; i < totalReadNumber; i++) {
	  bestHitMappingInfo[i].loc = -1;
	  bestHitMappingInfo[i].tprob = 0.0;
	}
}

void finalizeBestSingleMapping()
{
	int i = 0;
	char *_tmpQual, *_tmpSeq;
	char rqual[SEQ_LENGTH + 1];
	rqual[SEQ_LENGTH] = '\0';

	for(i = 0; i < _msf_seqListSize; i++) {
	  if(_msf_seqList[i].hits[0] != 0) {
		  if (bestHitMappingInfo[i].dir) {
			  reverse(_msf_seqList[i].qual, rqual, SEQ_LENGTH);
			  _tmpQual = rqual;
			  _tmpSeq = _msf_seqList[i].rseq;
		  }
		  else {
			  _tmpQual = _msf_seqList[i].qual;
			  _tmpSeq = _msf_seqList[i].seq;
		  }
	  _msf_output.QNAME = _msf_seqList[i].name;
	  _msf_output.FLAG = 16 * bestHitMappingInfo[i].dir;
	  _msf_output.RNAME = bestHitMappingInfo[i].chr;
	  _msf_output.POS = bestHitMappingInfo[i].loc;

	  _msf_output.MAPQ = mapQ(i);

	  _msf_output.CIGAR = bestHitMappingInfo[i].cigar;
	  _msf_output.MRNAME = "*";
	  _msf_output.MPOS = 0;
	  _msf_output.ISIZE = 0;

	  _msf_output.SEQ = _tmpSeq;
	  _msf_output.QUAL = _tmpQual;

	  _msf_output.optSize = 2;
	  _msf_output.optFields = _msf_optionalFields;

	  _msf_optionalFields[0].tag = "NM";
	  _msf_optionalFields[0].type = 'i';
	  _msf_optionalFields[0].iVal = bestHitMappingInfo[i].err;

	  _msf_optionalFields[1].tag = "MD";
	  _msf_optionalFields[1].type = 'Z';
	  _msf_optionalFields[1].sVal = bestHitMappingInfo[i].md;

	  output(_msf_output);
	}
	}
	freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}

int compare(const void *a, const void *b) {
	return ((Pair *) a)->hv - ((Pair *) b)->hv;
}

void preProcessReads() {
	int i = 0;
	//printf("reads=%d\n", _msf_seqListSize);
	_msf_sort_seqList = getMem(_msf_seqListSize * sizeof(Pair));
	for (i = 0; i < _msf_seqListSize; i++) {
		_msf_sort_seqList[i].hv = hashVal(_msf_seqList[i].seq);
		_msf_sort_seqList[i].readNumber = i;
	//printf("%d", _msf_sort_seqList[i].readNumber);
	}
	qsort(_msf_sort_seqList, _msf_seqListSize, sizeof(Pair), compare);
	_msf_map_sort_seqList = getMem(_msf_seqListSize * sizeof(int));
	for (i = 0; i < _msf_seqListSize; i++)
		_msf_map_sort_seqList[_msf_seqList[i].readNumber] = i;
}

void initFAST(Read *seqList, int seqListSize, int *samplingLocs,
	      int samplingLocsSize, char *genFileName) {
	  int i;
	  if (_msf_optionalFields == NULL) {
	      _msf_op = getMem(SEQ_LENGTH);
	      if (pairedEndMode) {
	          _msf_optionalFields = getMem(8 * sizeof(OPT_FIELDS));
	      } else {
	          _msf_optionalFields = getMem(2 * sizeof(OPT_FIELDS));
	      }

	      for (i = 0; i < 200; i++) {
	          sprintf(_msf_numbers[i], "%d%c", i, '\0');
	      }
	      sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
	  }

	  if (_msf_samplingLocsEnds == NULL) {
	      _msf_samplingLocs = samplingLocs;
	      _msf_samplingLocsSize = samplingLocsSize;

	      _msf_samplingLocsEnds = getMem(sizeof(int) * _msf_samplingLocsSize);
	      for (i = 0; i < _msf_samplingLocsSize; i++) {
	          _msf_samplingLocsEnds[i] = _msf_samplingLocs[i] + WINDOW_SIZE - 1;
	      }

	      _msf_seqList = seqList;
	      _msf_seqListSize = seqListSize;
	      preProcessReads();
	      _msf_oeaMapping = getMem(_msf_seqListSize * sizeof(int));
	      for (i = 0; i < _msf_seqListSize; i++) {
	          _msf_oeaMapping[i] = 0;
	      }

	      _msf_discordantMapping = getMem(_msf_seqListSize * sizeof(int));
	      for (i = 0; i < _msf_seqListSize; i++) {
	          _msf_discordantMapping[i] = 0;
	      }

	  }

	  if (_msf_refGenName == NULL) {
	      _msf_refGenName = getMem(4 * SEQ_LENGTH);
	  }
	  _msf_refGen = getRefGenome();
	  _msf_refGenLength = strlen(_msf_refGen);

	  _msf_refGenOffset = getRefGenomeOffset();
	  snprintf(_msf_refGenName, 4 * SEQ_LENGTH, "%s%c", getRefGenomeName(), '\0');
	  _msf_refGenName[strlen(getRefGenomeName())] = '\0';

	  if (_msf_verifiedLocs != NULL) {
	      freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength + 1));
	  }
	  _msf_verifiedLocs = (int *) getMem(sizeof(int) * (_msf_refGenLength + 1));

	  for (i = 0; i <= _msf_refGenLength; i++)
	      _msf_verifiedLocs[i] = _msf_seqListSize * 10 + 1;


	  if (pairedEndMode && _msf_seqHits == NULL) {
	      _msf_mappingInfo = getMem(seqListSize * sizeof(MappingInfo));
	      for (i = 0; i < seqListSize; i++) {
	          _msf_mappingInfo[i].next = NULL;
	          _msf_mappingInfo[i].size = 0;
	      }
	      _msf_seqHits = getMem((_msf_seqListSize) * sizeof(int));
	      for (i = 0; i < _msf_seqListSize; i++) {
	          _msf_seqHits[i] = 0;
	      }
	      _msf_readHasConcordantMapping = getMem(
	              _msf_seqListSize / 2 * sizeof(char));
	      for (i = 0; i < _msf_seqListSize / 2; i++) {
	          _msf_readHasConcordantMapping[i] = 0;
	      }
	      initLoadingRefGenome(genFileName);
	  }

	  if (_msf_refGenOffset == 0) {
	      _msf_refGenBeg = 1;
	  } else {
	      _msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2 + errThreshold;
	  }
	  _msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;
}

void finalizeFAST() {
	freeMem(_msf_seqHits, (_msf_seqListSize) * sizeof(int));
	freeMem(_msf_refGenName, 4 * SEQ_LENGTH);
	//freeMem(_msf_map_sort_seqList, sizeof(Pair) * _msf_seqListSize);
	//freeMem(_msf_sort_seqList, sizeof(int) * _msf_seqListSize);
}

int addCigarSize(int cnt) {
	if (cnt < 10)
	  return 1;
	else if (cnt < 100)
	  return 2;
	return 3;
}

// Generate Cigar from the back tracking matrix
void generateCigar(char *matrix, int matrixLength, char *cigar) {

	int i = 0;
	int counterM = 0;
	int counterI = 0;
	int counterD = 0;
	int cigarSize = 0;

	cigar[0] = '\0';
	while (i < matrixLength) {
		if (matrix[i] == 'M') {
			counterM++;
			if (counterI != 0) {
				sprintf(cigar, "%s%dI", cigar, counterI);
				cigarSize += addCigarSize(counterI) + 1;
				cigar[cigarSize] = '\0';
				counterI = 0;
			} else if (counterD != 0) {
				sprintf(cigar, "%s%dD", cigar, counterD);
				cigarSize += addCigarSize(counterD) + 1;
				cigar[cigarSize] = '\0';
				counterD = 0;
			}
		} else if (matrix[i] == 'I') {
			if (counterM != 0) {
				sprintf(cigar, "%s%dM", cigar, counterM);
				cigarSize += addCigarSize(counterM) + 1;
				cigar[cigarSize] = '\0';
				counterM = 0;
			} else if (counterD != 0) {
				sprintf(cigar, "%s%dD", cigar, counterD);
				cigarSize += addCigarSize(counterD) + 1;
				cigar[cigarSize] = '\0';
				counterD = 0;
			}
			counterI++;
			i++;
		} else if (matrix[i] == 'D') {
			if (counterM != 0) {
				sprintf(cigar, "%s%dM", cigar, counterM);
				cigarSize += addCigarSize(counterM) + 1;
				cigar[cigarSize] = '\0';
				counterM = 0;
			} else if (counterI != 0) {
				sprintf(cigar, "%s%dI", cigar, counterI);
				cigarSize += addCigarSize(counterI) + 1;
				cigar[cigarSize] = '\0';
				counterI = 0;
			}
			counterD++;
			i++;
		} else {
			counterM++;
			if (counterI != 0) {
				sprintf(cigar, "%s%dI", cigar, counterI);
				cigarSize += addCigarSize(counterI) + 1;
				cigar[cigarSize] = '\0';
				counterI = 0;
			} else if (counterD != 0) {
				sprintf(cigar, "%s%dD", cigar, counterD);
				cigarSize += addCigarSize(counterD) + 1;
				cigar[cigarSize] = '\0';
				counterD = 0;
			}
		}
		i++;
	}
	if (counterM != 0) {
		sprintf(cigar, "%s%dM", cigar, counterM);
		cigarSize += addCigarSize(counterM) + 1;
		cigar[cigarSize] = '\0';
		counterM = 0;
	} else if (counterI != 0) {
		sprintf(cigar, "%s%dI", cigar, counterI);
		cigarSize += addCigarSize(counterI) + 1;
		cigar[cigarSize] = '\0';
		counterI = 0;
	} else if (counterD != 0) {
		sprintf(cigar, "%s%dD", cigar, counterD);
		cigarSize += addCigarSize(counterD) + 1;
		cigar[cigarSize] = '\0';
		counterD = 0;
	}
	cigar[cigarSize] = '\0';
}

// Creates the Cigar output from the mismatching positions format  [0-9]+(([ACTGN]|\^[ACTGN]+)[0-9]+)*
void generateCigarFromMD(char *mismatch, int mismatchLength, char *cigar) {
	int i = 0;
	int j = 0;

	int start = 0;
	int cigarSize = 0;

	cigar[0] = '\0';

	while (i < mismatchLength) {
		if (mismatch[i] >= '0' && mismatch[i] <= '9') {
			start = i;
			while (mismatch[i] >= '0' && mismatch[i] <= '9'
					&& i < mismatchLength)
				i++;

			int value = atoi(mismatch + start);
			for (j = 0; j < value - 1; j++) {
				cigar[cigarSize] = 'M';
				cigarSize++;
			}
			cigar[cigarSize] = 'M';
		} else if (mismatch[i] == '^') {
			cigar[cigarSize] = 'I';
			i++;
		} else if (mismatch[i] == '\'') {
			cigar[cigarSize] = 'D';
			i++;
		} else {
			cigar[cigarSize] = 'M';
			cigarSize++;
		}
		cigarSize++;
		i++;
	}
	cigar[cigarSize] = '\0';
}

void generateSNPSAM(char *matrix, int matrixLength, char *outputSNP) {

	int i = 0;

	int counterM = 0;
	int counterD = 0;

	char delete[100];

	int snpSize = 0;

	outputSNP[0] = '\0';
	delete[0] = '\0';

	while (i < matrixLength) {
	  if (matrix[i] == 'M') {
	    counterM++;
	    if (counterD != 0) {
	delete[counterD] = '\0';
	counterD = 0;
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	delete[0] = '\0';
	    }
	  } else if (matrix[i] == 'D') {
	    if (counterM != 0) {
	sprintf(outputSNP, "%s%d", outputSNP, counterM);
	snpSize += addCigarSize(counterM);
	outputSNP[snpSize] = '\0';
	counterM = 0;
	delete[counterD] = matrix[i + 1];
	i++;
	counterD++;
	    } else if (counterD != 0) {
	delete[counterD] = matrix[i + 1];
	counterD++;
	i++;
	    } else {
	delete[counterD] = matrix[i + 1];
	counterD++;
	i++;
	    }
	  } else if (matrix[i] == 'I') {
	    if (counterM != 0) {
	// sprintf(outputSNP, "%s%d\0", outputSNP, counterM);
	//counterM++;
	    } else if (counterD != 0) {
	delete[counterD] = '\0';
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	counterD = 0;
	delete[0] = '\0';
	    }
	    i++;

	  } else {
	    if (counterM != 0) {
	sprintf(outputSNP, "%s%d", outputSNP, counterM);
	snpSize += addCigarSize(counterM);
	outputSNP[snpSize] = '\0';
	counterM = 0;
	    }
	    if (counterD != 0) {
	delete[counterD] = '\0';
	counterD = 0;
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	delete[0] = '\0';
	    }
	    sprintf(outputSNP, "%s%c", outputSNP, matrix[i]);
	    snpSize += 1;
	    outputSNP[snpSize] = '\0';
	  }
	  i++;
	}

	if (counterM != 0) {
	  sprintf(outputSNP, "%s%d", outputSNP, counterM);
	  snpSize += addCigarSize(counterM);
	  outputSNP[snpSize] = '\0';
	  counterM = 0;
	} else if (counterD != 0) {
	  delete[counterD] = '\0';
	  sprintf(outputSNP, "%s^%s", outputSNP, delete);
	  snpSize += strlen(delete) + 1;
	  outputSNP[snpSize] = '\0';
	  counterD = 0;
	}

	outputSNP[snpSize] = '\0';
}

int compareEntrySize(const void *a, const void *b) {
	return ((*(key_struct*) a).key_entry_size
	  - (*(key_struct*) b).key_entry_size);
}

int compareOut(const void *a, const void *b) {
	FullMappingInfo *aInfo = (FullMappingInfo *) a;
	FullMappingInfo *bInfo = (FullMappingInfo *) b;
	return aInfo->loc - bInfo->loc;
}

void outputPairFullMappingInfo(FILE *fp, int readNumber) {

	char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
	char rqual1[SEQ_LENGTH + 1], rqual2[SEQ_LENGTH + 1];

	rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';

	seq1 = _msf_seqList[readNumber * 2].seq;
	rseq1 = _msf_seqList[readNumber * 2].rseq;
	qual1 = _msf_seqList[readNumber * 2].qual;

	reverse(_msf_seqList[readNumber * 2].qual, rqual1, SEQ_LENGTH);

	seq2 = _msf_seqList[readNumber * 2 + 1].seq;
	rseq2 = _msf_seqList[readNumber * 2 + 1].rseq;
	qual2 = _msf_seqList[readNumber * 2 + 1].qual;

	reverse(_msf_seqList[readNumber * 2 + 1].qual, rqual2, SEQ_LENGTH);

	if (bestHitMappingInfo[readNumber * 2].loc == -1
	    && bestHitMappingInfo[readNumber * 2 + 1].loc == -1)
	  return;
	else {
	  char *seq;
	  char *qual;
	  char d1;
	  char d2;
	  int isize;
	  int proper = 0;
	  // ISIZE CALCULATION
	  // The distance between outer edges
	  isize = abs(
		bestHitMappingInfo[readNumber * 2].loc
		- bestHitMappingInfo[readNumber * 2 + 1].loc)
	    + SEQ_LENGTH - 2;

	  if (bestHitMappingInfo[readNumber * 2].loc
	- bestHitMappingInfo[readNumber * 2 + 1].loc > 0) {
	    isize *= -1;
	  }
	  d1 = (bestHitMappingInfo[readNumber * 2].dir == -1) ? 1 : 0;
	  d2 = (bestHitMappingInfo[readNumber * 2 + 1].dir == -1) ? 1 : 0;
	  if (d1) {
	    seq = rseq1;
	    qual = rqual1;
	  } else {
	    seq = seq1;
	    qual = qual1;
	  }
	  //TODO for CG like SOLID
	  if ( (d1 && d2) || (!d1 && !d2)) {
	    proper = 2;
	  } else {
	    proper = 0;
	  }

	  _msf_output.POS = bestHitMappingInfo[readNumber * 2].loc;
	  _msf_output.MPOS = bestHitMappingInfo[readNumber * 2 + 1].loc;
	  _msf_output.FLAG = 1 + proper + 16 * d1 + 32 * d2 + 64;
	  _msf_output.ISIZE = isize;
	  _msf_output.SEQ = seq;
	  _msf_output.QUAL = qual;
	  _msf_output.QNAME = _msf_seqList[readNumber * 2].name;
	  _msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;

	  _msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);

	  _msf_output.CIGAR = bestHitMappingInfo[readNumber * 2].cigar;
	  _msf_output.MRNAME = "=";

	  _msf_output.optSize = 2;
	  _msf_output.optFields = _msf_optionalFields;

	  _msf_optionalFields[0].tag = "NM";
	  _msf_optionalFields[0].type = 'i';
	  _msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber * 2].err;

	  _msf_optionalFields[1].tag = "MD";
	  _msf_optionalFields[1].type = 'Z';
	  _msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2].md;


	  output(_msf_output);
	  if (d2) {
	    seq = rseq2;
	    qual = rqual2;
	  } else {
	    seq = seq2;
	    qual = qual2;
	  }
	  _msf_output.POS = bestHitMappingInfo[readNumber * 2 + 1].loc;
	  _msf_output.MPOS = bestHitMappingInfo[readNumber * 2].loc;
	  _msf_output.FLAG = 1 + proper + 16 * d2 + 32 * d1 + 128;
	  _msf_output.ISIZE = -isize;
	  _msf_output.SEQ = seq;
	  _msf_output.QUAL = qual;
	  _msf_output.QNAME = _msf_seqList[readNumber * 2].name;
	  _msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;

	  _msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);

	  _msf_output.CIGAR = bestHitMappingInfo[readNumber * 2 + 1].cigar;
	  _msf_output.MRNAME = "=";

	  _msf_output.optSize = 2;
	  _msf_output.optFields = _msf_optionalFields;

	  _msf_optionalFields[0].tag = "NM";
	  _msf_optionalFields[0].type = 'i';
	  _msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber * 2 + 1].err;

	  _msf_optionalFields[1].tag = "MD";
	  _msf_optionalFields[1].type = 'Z';
	  _msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2 + 1].md;
	  output(_msf_output);
	}
}

/*
	Find the closet one to the c
	@return 0: if the x1 is closer to c
	1: if the x2 is closer to c
	2: if both distance are equal
	-1: if error
*/

int findNearest(int x1, int x2, int c) {
	if (abs(x1 - c) < abs(x2 - c))
	  return 0;
	else if (abs(x1 - c) > abs(x2 - c))
	  return 1;
	else if (abs(x1 - c) == abs(x2 - c))
	  return 2;
	else
	  return -1;
}

//void initBestMapping(int totalReadNumber) {
//	bestHitMappingInfo = getMem(totalReadNumber * sizeof(BestFullMappingInfo));
//
//	int i = 0;
//	for (i = 0; i < totalReadNumber; i++) {
//	  bestHitMappingInfo[i].loc = -1;
//	  bestHitMappingInfo[i].tprob = 0.0;
//	}
//}

void performBestMapping(int readCnt) {
	int i = 0;
	for (i = 0; i < readCnt / 2; i++)
		outputPairFullMappingInfo(NULL, i);
}

void finalizeBestMapping(int readCnt) {
	freeMem(bestHitMappingInfo, readCnt * sizeof(BestFullMappingInfo));
}

void finalizeBestConcordantDiscordant() {
	  int i = 0;
	  for (i = 0; i < _msf_seqListSize / 2; i++) {
	      outputPairFullMappingInfo(NULL, i);
	  }
	  freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(BestFullMappingInfo));
}

double mapProb(int readNumber, char *md, int dir, int err) {
	int i = 0;
	int mdlen = strlen(md);
	char buf[MAX_CIGAR_SIZE];
	int j = 0;
	double phred = 0.0;
	int errloc = 0;
	int errcnt = 0; //since I cannot calculate deletion base quality
	buf[0] = 0;

	if (err == 0)
	  return 1.0;

	while (i<mdlen){
	  if (isdigit(md[i])) buf[j++]=md[i++];
	  else if (isalpha(md[i])){
	    /* mismatch */
	    errcnt++;
	    buf[j] = '\0';
	    if (j != 0)
	errloc += atoi(buf);
	    else if (i!=0)
	errloc++;

	    j=0; buf[0]=0;

	    if (dir)
	phred += (double) (_msf_seqList[readNumber].qual[SEQ_LENGTH-errloc-1] - 33);
	    else
	phred += (double) (_msf_seqList[readNumber].qual[errloc] - 33);

	    i++;
	  }

	  else if (md[i]=='^'){
	    /* insertion to the read / deletion from reference  */
	    if (j!=0){
	buf[j]=0;
	errloc += atoi(buf);
	buf[0] = 0;
	    }
	    j=0;
	    i++; /* pass ^ */
	    while (isalpha(md[i++])) j++;
	    errloc += j;
	    j = 0;
	  }
	}

	double indel_prob = 1;
	if (errcnt != err)
	  indel_prob = 0.0002 * (err - errcnt);

	return pow(10, -1 * (phred / 10)) * indel_prob;

}

int mapQ(int readNumber) {
	  int mapqual;
	  double mapprob;
	  mapprob = mapProb(readNumber, bestHitMappingInfo[readNumber].md,
	          bestHitMappingInfo[readNumber].dir, bestHitMappingInfo[readNumber].err);
	  if (mapprob == bestHitMappingInfo[readNumber].tprob)
	      mapqual = 40;
	  else
	      mapqual =  (int) (round(-10.0 * log10(1 - (mapprob / bestHitMappingInfo[readNumber].tprob))));
	  if (mapqual > 40) mapqual = 40;
	  return mapqual;
}

void setFullMappingInfo(int readNumber, int loc, int dir, int err, int score,
			char *md, char * refName, char *cigar) {
	bestHitMappingInfo[readNumber].loc = loc;
	bestHitMappingInfo[readNumber].dir = dir;
	bestHitMappingInfo[readNumber].err = err;
	bestHitMappingInfo[readNumber].score = score;

	strncpy(bestHitMappingInfo[readNumber].md, md, strlen(md) + 1);
	strncpy(bestHitMappingInfo[readNumber].chr, refName, strlen(refName) + 1);
	strncpy(bestHitMappingInfo[readNumber].cigar, cigar, strlen(cigar) + 1);
}

void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1,
			    FullMappingInfo mi2) {

	bestHitMappingInfo[readNumber * 2].loc = mi1.loc;
	bestHitMappingInfo[readNumber * 2].dir = mi1.dir;
	bestHitMappingInfo[readNumber * 2].err = mi1.err;
	bestHitMappingInfo[readNumber * 2].score = mi1.score;
	snprintf(bestHitMappingInfo[readNumber * 2].chr, MAX_REF_SIZE, "%s",
	   _msf_refGenName);

	strncpy(bestHitMappingInfo[readNumber * 2].md, mi1.md, strlen(mi1.md) + 1);
	strncpy(bestHitMappingInfo[readNumber * 2].cigar, mi1.cigar,
	  strlen(mi1.cigar) + 1);

	bestHitMappingInfo[readNumber * 2 + 1].loc = mi2.loc;
	bestHitMappingInfo[readNumber * 2 + 1].dir = mi2.dir;
	bestHitMappingInfo[readNumber * 2 + 1].err = mi2.err;
	bestHitMappingInfo[readNumber * 2 + 1].score = mi2.score;

	snprintf(bestHitMappingInfo[readNumber * 2 + 1].chr, MAX_REF_SIZE, "%s",
	   _msf_refGenName);

	strncpy(bestHitMappingInfo[readNumber * 2 + 1].md, mi2.md,
	  strlen(mi2.md) + 1);
	strncpy(bestHitMappingInfo[readNumber * 2 + 1].cigar, mi2.cigar,
	  strlen(mi2.cigar) + 1);

}

int outputPairedEnd(int pre_unmappedCnt) {
	int i = 0;
	char cigar[MAX_CIGAR_SIZE];
	int tmpOut;
	FILE* in1[_msf_openFiles];
	FILE* in2[_msf_openFiles];
	char fname1[_msf_openFiles][FILE_NAME_LENGTH];
	char fname2[_msf_openFiles][FILE_NAME_LENGTH];

	FILE *out = NULL, *out1 = NULL;
	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	int meanDistanceMapping = 0;
	char rqual1[SEQ_LENGTH + 1];
	char rqual2[SEQ_LENGTH + 1];
	int tmp = 0;

	loadRefGenome(&_msf_refGen, &_msf_refGenName, &tmpOut);

	if (pairedEndDiscordantMode) {
	  sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
	  sprintf(fname4, "%s__%s__oea", mappingOutputPath, mappingOutput);
	  out = fileOpen(fname3, "a");
	  out1 = fileOpen(fname4, "a");
	}
	FullMappingInfo *mi1 = getMem(sizeof(FullMappingInfo) * _msf_maxLSize);
	FullMappingInfo *mi2 = getMem(sizeof(FullMappingInfo) * _msf_maxRSize);
	_msf_fileCount[_msf_maxFile] = 0;
	for (i = 0; i < _msf_openFiles; i++) {
	  sprintf(fname1[i], "%s__%s__%s__%d__1.tmp", mappingOutputPath,
	    _msf_refGenName, mappingOutput, i);
	  sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][0],
	    "%s", fname1[i]);
	  sprintf(fname2[i], "%s__%s__%s__%d__2.tmp", mappingOutputPath,
	    _msf_refGenName, mappingOutput, i);
	  sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][1],
	    "%s", fname2[i]);
	  in1[i] = fileOpen(fname1[i], "r");
	  in2[i] = fileOpen(fname2[i], "r");
	  _msf_fileCount[_msf_maxFile]++;
	}
	_msf_maxFile++;
	int size;
	int j, k;
	int size1, size2;
	meanDistanceMapping =
	  (pairedEndDiscordantMode == 1) ?
	  (minPairEndedDiscordantDistance
	   + maxPairEndedDiscordantDistance) / 2 + SEQ_LENGTH :
	  (minPairEndedDistance + maxPairEndedDistance) / 2
	  + SEQ_LENGTH;

	for (i = 0; i < _msf_seqListSize / 2; i++) {
	    size1 = size2 = 0;

	    for (j = 0; j < _msf_openFiles; j++) {
	        tmpOut = fread(&size, sizeof(int), 1, in1[j]);
	        if (size > 0) {
	            for (k = 0; k < size; k++) {
	                mi1[size1 + k].dir = 1;
	                tmpOut = fread(&(mi1[size1 + k].loc), sizeof(int), 1, in1[j]);
	                tmpOut = fread(&(mi1[size1 + k].err), sizeof(int), 1, in1[j]);
	                tmpOut = fread(&(mi1[size1 + k].cigarSize), sizeof(int), 1, in1[j]);
	                tmpOut = fread((mi1[size1 + k].cigar), sizeof(char), mi1[size1 + k].cigarSize, in1[j]);
	                mi1[size1 + k].cigar[mi1[size1 + k].cigarSize] = '\0';
	                tmpOut = fread(&(mi1[size1 + k].mdSize), sizeof(int), 1, in1[j]);
	                tmpOut = fread((mi1[size1 + k].md), sizeof(char), (mi1[size1 + k].mdSize), in1[j]);
	                mi1[size1 + k].md[mi1[size1 + k].mdSize] = '\0';

	                if (mi1[size1 + k].loc < 1) {
	                    mi1[size1 + k].loc *= -1;
	                    mi1[size1 + k].dir = -1;
	                }
	            }
	            qsort(mi1 + size1, size, sizeof(FullMappingInfo), compareOut);
	            size1 += size;
	        }
	    }

	    for (j = 0; j < _msf_openFiles; j++) {
	        tmpOut = fread(&size, sizeof(int), 1, in2[j]);
	        if (size > 0) {
	            for (k = 0; k < size; k++) {
	                mi2[size2 + k].dir = 1;
	                tmpOut = fread(&(mi2[size2 + k].loc), sizeof(int), 1, in2[j]);
	                tmpOut = fread(&(mi2[size2 + k].err), sizeof(int), 1, in2[j]);
	                tmpOut = fread(&(mi2[size2 + k].cigarSize), sizeof(int), 1, in2[j]);
	                tmpOut = fread((mi2[size2 + k].cigar), sizeof(char), mi2[size2 + k].cigarSize, in2[j]);
	                mi2[size2 + k].cigar[mi2[size2 + k].cigarSize] = '\0';
	                tmpOut = fread(&(mi2[size2 + k].mdSize), sizeof(int), 1, in2[j]);
	                tmpOut = fread((mi2[size2 + k].md), sizeof(char), mi2[size2 + k].mdSize, in2[j]);
	                mi2[size2 + k].md[mi2[size2 + k].mdSize] = '\0';
	                if (mi2[size2 + k].loc < 1) {
	                    mi2[size2 + k].loc *= -1;
	                    mi2[size2 + k].dir = -1;
	                }
	            }
	            qsort(mi2 + size2, size, sizeof(FullMappingInfo), compareOut);
	            size2 += size;
	        }
	    }


	  int lm, ll, rl, rm;
	  int pos = 0;

	  if (pairedEndDiscordantMode) {
	      for (j = 0; j < size1; j++) {
	          lm = mi1[j].loc - maxPairEndedDiscordantDistance + 1;
	          ll = mi1[j].loc - minPairEndedDiscordantDistance + 1;
	          rl = mi1[j].loc + minPairEndedDiscordantDistance - 1;
	          rm = mi1[j].loc + maxPairEndedDiscordantDistance - 1;
	          while (pos < size2 && mi2[pos].loc < lm) {
	              pos++;
	          }
	          k = pos;
	          while (k < size2 && mi2[k].loc <= rm) {
	              if (mi2[k].loc <= ll || mi2[k].loc >= rl) {
	                  if (  (mi1[j].dir == 1 && mi2[k].dir == 1)
	                          || (mi1[j].dir == -1 && mi2[k].dir == -1)) {
	                      _msf_seqList[i * 2].hits[0] = 1;
	                      _msf_seqList[i * 2 + 1].hits[0] = 1;
	                      if (nosamMode != 0) {
	                          size1 = 0;
	                          size2 = 0;
	                      }
	                      break;
	                  }
	              }
	              k++;
	          }
	      }
	      _msf_seqHits[i * 2] += size1;
	      _msf_seqHits[i * 2 + 1] += size2;
	      if (_msf_seqHits[i * 2 + 1] * _msf_seqHits[i * 2]
	              > DISCORDANT_CUT_OFF && nosamMode != 0) {
	          _msf_seqList[i * 2].hits[0] = 1;
	          _msf_seqList[i * 2 + 1].hits[0] = 1;
	          size1 = 0;
	          size2 = 0;
	      }
	    int rNo = 0;
	    int loc = 0;
	    int err = 0;
	    float sc = 0;
	    char l = 0;
	    //write the OEA data
	    if (_msf_seqHits[i * 2] == 0){
	        for (k = 0; k < size2 && _msf_oeaMapping[i * 2 + 1] < maxOEAOutput; k++) {
	            rNo = i * 2 + 1;
	            loc = mi2[k].loc * mi2[k].dir;
	            err = mi2[k].err;
	            sc = mi2[k].score;
	            l = strlen(_msf_refGenName);
	            tmp = fwrite(&rNo, sizeof(int), 1, out1);
	            tmp = fwrite(&l, sizeof(char), 1, out1);
	            tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);
	            tmp = fwrite(&loc, sizeof(int), 1, out1);
	            tmp = fwrite(&err, sizeof(int), 1, out1);
	            tmp = fwrite(&sc, sizeof(float), 1, out1);
	            if (mi2[k].cigarSize > SEQ_LENGTH || mi2[k].cigarSize <= 0)
	                printf("ERROR  CIGAR size=%d %s\n", mi2[k].cigarSize, _msf_seqList[i * 2 + 1].seq);
	            tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1, out1);
	            tmp = fwrite((mi2[k].cigar), sizeof(char), mi2[k].cigarSize, out1);
	            tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out1);
	            //if (mi2[k].mdSize > SEQ_LENGTH || mi1[k].mdSize <= 0)
	            if (mi2[k].mdSize > SEQ_LENGTH || mi1[k].mdSize < 0)
	                printf("ERROR mdsize = %d %s\n", mi1[k].mdSize, _msf_seqList[i * 2].seq);
	            tmp = fwrite((mi2[k].md), sizeof(char), mi2[k].mdSize, out1);
	            _msf_oeaMapping[i * 2 + 1]++;
	        }
	    }
	    if (_msf_seqHits[i * 2 + 1] == 0){
	        for (j = 0; j < size1 && _msf_oeaMapping[i * 2] < maxOEAOutput; j++) {
	            rNo = i * 2;
	            loc = mi1[j].loc * mi1[j].dir;
	            err = mi1[j].err;
	            sc = mi1[j].score;
	            l = strlen(_msf_refGenName);
	            tmp = fwrite(&rNo, sizeof(int), 1, out1);
	            tmp = fwrite(&l, sizeof(char), 1, out1);
	            tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);
	            tmp = fwrite(&loc, sizeof(int), 1, out1);
	            tmp = fwrite(&err, sizeof(int), 1, out1);
	            tmp = fwrite(&sc, sizeof(float), 1, out1);
	            if (mi1[j].cigarSize > SEQ_LENGTH || mi1[j].cigarSize <= 0)
	                printf("ERROR %d %s\n", mi1[j].cigarSize, _msf_seqList[i * 2 + 1].seq);
	            tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1, out1);
	            tmp = fwrite((mi1[j].cigar), sizeof(char), mi1[j].cigarSize, out1);
	            tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out1);
	            //if (mi1[j].mdSize > SEQ_LENGTH || mi1[j].mdSize <= 0)
	            if (mi1[j].mdSize > SEQ_LENGTH || mi1[j].mdSize < 0)
	                printf("ERROR mdsize = %d %s\n", mi1[j].mdSize, _msf_seqList[i * 2 + 1].seq);
	            tmp = fwrite((mi1[j].md), sizeof(char), mi1[j].mdSize, out1);
	            _msf_oeaMapping[i * 2]++;
	        }
	    }
	  }
	  char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
	  rqual1[SEQ_LENGTH] = '\0';
	  rqual2[SEQ_LENGTH] = '\0';
	  rqual1[0] = '\0';
	  rqual2[0] = '\0';

	  seq1 = _msf_seqList[i * 2].seq;
	  rseq1 = _msf_seqList[i * 2].rseq;
	  qual1 = _msf_seqList[i * 2].qual;
	  strncpy(rqual1, _msf_seqList[i * 2].qual, SEQ_LENGTH);
	  seq2 = _msf_seqList[i * 2 + 1].seq;
	  rseq2 = _msf_seqList[i * 2 + 1].rseq;
	  qual2 = _msf_seqList[i * 2 + 1].qual;
	  strncpy(rqual2, _msf_seqList[i * 2 + 1].qual, SEQ_LENGTH);
	  if (pairedEndDiscordantMode) {
	      for (k = 0; k < size1; k++) {
	          mi1[k].score = calculateScore(mi1[k].loc,
	                  (mi1[k].dir == -1) ? rseq1 : seq1,
	                  (mi1[k].dir == -1) ? rqual1 : qual1, mi1[k].cigar);
	      }
	      for (k = 0; k < size2; k++) {
	          mi2[k].score = calculateScore(mi2[k].loc,
	                  (mi2[k].dir == -1) ? rseq2 : seq2,
	                  (mi2[k].dir == -1) ? rqual2 : qual2, mi2[k].cigar);
	      }
	  }
	  // CALKAN MAPQ FOR PE
	  for (j = 0; j < size1; j++) {
	          if (mi1[j].err != 0){
	              bestHitMappingInfo[i*2].tprob += mapProb(i*2, mi1[j].md, mi1[j].dir, mi1[j].err);
	          }
	      }
	      for (k = 0; k < size2; k++) {
	          if (mi2[k].err != 0){
	      bestHitMappingInfo[i*2+1].tprob += mapProb((i*2+1), mi2[k].md, mi2[k].dir, mi2[k].err);
	          }
	      }

		int total_num = 0;	// DHL Temp
	  if (pairedEndDiscordantMode) {
	      for (j = 0; j < size1; j++) {
	          for (k = 0; k < size2; k++) {
	              int dir1 = mi1[j].dir;
	              int dir2 = mi2[k].dir;
	              int loc1 = mi1[j].loc;
	              int loc2 = mi2[k].loc;
	              int best_err1 = bestHitMappingInfo[i * 2].err;
	              int best_err2 = bestHitMappingInfo[i * 2+1].err;
	              int best_loc1 = bestHitMappingInfo[i * 2].loc;
	              int best_loc2 = bestHitMappingInfo[i * 2+1].loc;
	              if ( (( dir1 > 0 && dir2 > 0 ) || (dir1 < 0 && dir2 < 0))  &&
	                  (( loc1 != -1 || loc2 != -1)  &&
	                  (abs(loc1 - loc2) >  minPairEndedDiscordantDistance) &&
	                  (abs(loc1 - loc2) <  maxPairEndedDiscordantDistance))) {
	                  //POSSIBLE CONCORDANT
	                  if(_msf_readHasConcordantMapping[i] == 0) {
	                      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                      _msf_readHasConcordantMapping[i] = 1;
	                      _msf_seqList[i * 2].hits[0] = 1;
	                      _msf_seqList[i * 2 + 1].hits[0] = 1;
	                  } else {
	                      if (best_err1+best_err2 >= mi1[j].err + mi2[k].err) {
	                          if ( best_err1+best_err2 == mi1[j].err + mi2[k].err
	                             && findNearest(
	                                      abs(best_loc1-best_loc2),
	                                      abs(loc2-loc1),
	                                      meanDistanceMapping) == 0) {
	                              continue;
	                          }
	                          setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                      }
	                  }
	              }
	              //DISCORDANT TO TEMP FILE FOR POST PROCESSING
	              else if (_msf_readHasConcordantMapping[i] == 0
	                      && _msf_seqHits[i * 2] != 0
	                      && _msf_seqHits[i * 2 + 1] != 0) {
	                  int rNo = i;
	                  int loc = mi1[j].loc * mi1[j].dir;
	                  int err = mi1[j].err;
	                  float sc = mi1[j].score;
	                  char l = strlen(_msf_refGenName);
	                  if (_msf_discordantMapping[i * 2] < maxDiscordantOutput) {
	                      tmp = fwrite(&rNo, sizeof(int), 1, out);
	                      tmp = fwrite(&l, sizeof(char), 1, out);
	                      tmp = fwrite(_msf_refGenName, sizeof(char), l, out);
	                      tmp = fwrite(&loc, sizeof(int), 1, out);
	                      tmp = fwrite(&err, sizeof(int), 1, out);
	                      tmp = fwrite(&sc, sizeof(float), 1, out);
	                      tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1, out);
	                      tmp = fwrite((mi1[j].cigar), sizeof(char), mi1[j].cigarSize, out);
	                      tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out);
	                      tmp = fwrite((mi1[j].md), sizeof(char), mi1[j].mdSize, out);
	                      loc = mi2[k].loc * mi2[k].dir;
	                      err = mi2[k].err;
	                      sc = mi2[k].score;
	                      tmp = fwrite(&loc, sizeof(int), 1, out);
	                      tmp = fwrite(&err, sizeof(int), 1, out);
	                      tmp = fwrite(&sc, sizeof(float), 1, out);
	                      tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1, out);
	                      tmp = fwrite((mi2[k].cigar), sizeof(char), mi2[k].cigarSize, out);
	                      tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out);
	                      tmp = fwrite((mi2[k].md), sizeof(char), mi2[k].mdSize, out);
	                      _msf_discordantMapping[i * 2]++;
										total_num++;
	                  }
	                  //SET THE BEST DISCORDANT
	                  //BEGIN {Farhad Hormozdiari}
	                  if (best_loc1 == -1 && best_loc2 == -1 && _msf_readHasConcordantMapping[i] == 0) {
	                      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                      _msf_seqList[i * 2].hits[0] = 1;
	                      _msf_seqList[i * 2 + 1].hits[0] = 1;
	                  } else if (best_err1 + best_err2 >= mi1[j].err + mi2[k].err && _msf_readHasConcordantMapping[i] == 0) {
	                      if (best_err1 + best_err2 == mi1[j].err + mi2[k].err
	                              && findNearest(
	                                  abs(best_loc2-best_loc1),
	                                  abs(loc1 - loc2),
	                                  meanDistanceMapping) == 0) {
	                          continue;
	                      }
	                      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                  }
	                  //END {Farhad Hormozdiari}
	              }
	          }
	      }
	  } else {
	      for (j = 0; j < size1; j++) {
	          for (k = 0; k < size2; k++) {
	               int dir1 = mi1[j].dir;
	               int dir2 = mi2[k].dir;
	               int loc1 = mi1[j].loc;
	               int loc2 = mi2[k].loc;
	               int best_err1 = bestHitMappingInfo[i * 2].err;
	               int best_err2 = bestHitMappingInfo[i * 2+1].err;
	               int best_loc1 = bestHitMappingInfo[i * 2].loc;
	               int best_loc2 = bestHitMappingInfo[i * 2+1].loc;

	               if ( abs (mi2[k].loc - mi1[j].loc) >= minPairEndedDistance
	                          && (abs(mi2[k].loc - mi1[j].loc) <= maxPairEndedDistance)
	                          && ( (dir1>0 && dir2>0) || ( dir1< 0 && dir2 <0) ) ) {
	                  char *seq;
	                  char *qual;
	                  char d1;
	                  char d2;
	                  int isize;
	                  int proper = 0;
	                  // ISIZE CALCULATION
	                  // The distance between outer edges
	                  isize = abs(mi1[j].loc - mi2[k].loc) + SEQ_LENGTH - 2;
	                  if (mi1[j].loc - mi2[k].loc > 0) {
	                      isize *= -1;
	                  }
	                  d1 = (mi1[j].dir == -1) ? 1 : 0;
	                  d2 = (mi2[k].dir == -1) ? 1 : 0;
	                  //SET THE READ HAS CONCORDANT MAPPING
	                  _msf_readHasConcordantMapping[i] = 1;
	                  if (d1) {
	                      seq = rseq1;
	                      qual = rqual1;
	                  } else {
	                      seq = seq1;
	                      qual = qual1;
	                  }
	                  if ((d1 && d2) || (!d1 && !d2)) {
	                      proper = 2;
	                  } else {
	                      proper = 0;
	                  }
	                  _msf_output.POS = mi1[j].loc;
	                  _msf_output.MPOS = mi2[k].loc;
	                  _msf_output.FLAG = 1 + proper + 16 * d1 + 32 * d2 + 64;
	                  _msf_output.ISIZE = isize;
	                  _msf_output.SEQ			= seq;
	                  _msf_output.QUAL		= qual;
	                  _msf_output.QNAME = _msf_seqList[i * 2].name;
	                  _msf_output.RNAME = _msf_refGenName;
	                  _msf_output.MAPQ = 255;
	                  _msf_output.CIGAR = cigar;
	                  _msf_output.MRNAME = "=";
	                  _msf_output.optSize = 2;
	                  _msf_output.optFields = _msf_optionalFields;
	                  _msf_optionalFields[0].tag = "NM";
	                  _msf_optionalFields[0].type = 'i';
	                  _msf_optionalFields[0].iVal = mi1[j].err;
	                  _msf_optionalFields[1].tag = "MD";
	                  _msf_optionalFields[1].type = 'Z';
	                  _msf_optionalFields[1].sVal = mi1[j].md;
	                  if (!bestMode)
	                      output(_msf_output);
	                  if (d2) {
	                      seq = rseq2;
	                      qual = rqual2;
	                  } else {
	                      seq = seq2;
	                      qual = qual2;
	                  }
	                  _msf_output.POS = mi2[k].loc;
	                  _msf_output.MPOS = mi1[j].loc;
	                  _msf_output.FLAG = 1 + proper + 16 * d2 + 32 * d1 + 128;
	                  _msf_output.ISIZE = -isize;
	                  _msf_output.SEQ		= seq;
	                  _msf_output.QUAL		= qual;
	                  _msf_output.QNAME = _msf_seqList[i * 2].name;
	                  _msf_output.RNAME = _msf_refGenName;
	                  _msf_output.MAPQ = 255;
	                  _msf_output.CIGAR = cigar;
	                  _msf_output.MRNAME = "=";
	                  _msf_output.optSize = 2;
	                  _msf_output.optFields = _msf_optionalFields;
	                  _msf_optionalFields[0].tag = "NM";
	                  _msf_optionalFields[0].type = 'i';
	                  _msf_optionalFields[0].iVal = mi2[k].err;
	                  _msf_optionalFields[1].tag = "MD";
	                  _msf_optionalFields[1].type = 'Z';
	                  _msf_optionalFields[1].sVal = mi2[k].md;
	                  if (!bestMode)
	                      output(_msf_output);
	                  //SET THE BEST CONCORDANT
	                  //BEGIN {Farhad Hormozdiari}
	                  if (best_loc1 == -1 && best_loc2 == -1) {
	                      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                  } else {
	                      if (best_err1 + best_err2 >= mi1[j].err + mi2[k].err) {
	                          if (best_err1+best_err2 == mi1[j].err + mi2[k].err
	                                  && findNearest(
	                                      abs(best_loc2 - best_loc1),
	                                      abs(loc2 - loc1),
	                                      meanDistanceMapping) == 0) {
	                              continue;
	                          }
	                          setPairFullMappingInfo(i, mi1[j], mi2[k]);
	                      }
	                  }
	                  //END   {Farhad Hormozdiari}
	              }
	          }
	      }
	  }
	}

	// ---- Verify End


	if (pairedEndDiscordantMode) {
		fclose(out);
		fclose(out1);
	}

	for (i = 0; i < _msf_openFiles; i++) {
		fclose(in1[i]);
		fclose(in2[i]);
		unlink(fname1[i]);
		unlink(fname2[i]);
	}
	tmp++;

	freeMem(mi1, sizeof(FullMappingInfo) * _msf_maxLSize);
	freeMem(mi2, sizeof(FullMappingInfo) * _msf_maxRSize);
	_msf_openFiles = 0;

	/* calkan counter */
	int unmappedCnt = 0;
	for (i = 0; i < _msf_seqListSize; i++) {
	  if (_msf_seqHits[i] == 0) unmappedCnt++;
	}
	mappedSeqCnt = _msf_seqListSize - unmappedCnt;

	// READ SLICE
	unmappedCnt = unmappedCnt + pre_unmappedCnt;
	mappedSeqCnt = _msf_totalSeqListSize - unmappedCnt;

	return unmappedCnt;
}

float str2int(char *str, int index1, int index2) {
	char tmp[SEQ_MAX_LENGTH];
	strncpy(tmp, &str[index1], index2 - index1);
	tmp[index2 - index1] = '\0';
	return atol(tmp);
}

double binomial_coefficient(int n, int k){
	double ret = 1.0;
	int i = 0;
	for (i=0; i<k; i++){
	  ret *= (n - i);
	  ret /= (k - i);
	}
	return ret;
}

float calculateScore(int index, char *seq, char *qual, char *md) {
	int i = 0;
	int j;
	char *ref;
	char *ver;
	float score = 1;
	char tmp[2 * SEQ_MAX_LENGTH];
	int value = 0;
	int end = 0;
	int index1 = 0;
	int index2 = 0;

	ref = _msf_refGen + index - 1;
	ver = seq;

	while (1) {
	    if (i >= strlen(md))
	        break;
	    index1 = i;
	    while (md[i] >= '0' && md[i] <= '9') {
	        i++;
	    }
	    index2 = i;
	    value = str2int(md, index1, index2);
	    if (md[i] == 'M') {
	        for (j = 0; j < value; j++) {
	            tmp[end] = 'M';
	            end++;
	        }
	    } else if (md[i] == 'I') {
	        for (j = 0; j < value; j++) {
	            tmp[end] = 'I';
	            end++;
	        }
	    } else if (md[i] == 'D') {
	        for (j = 0; j < value; j++) {
	            tmp[end] = 'D';
	            end++;
	        }
	    }
	    i++;
	}
	tmp[end] = '\0';
	j = 0;
	for (i = 0; i < end; i++) {
	    if (tmp[i] == 'M') {
	        if (*ref != *ver) {
	            score *= 0.001 + 1 / pow(10, ((qual[j] - 33) / 10.0));
	        }
	        ref++;
	        ver++;
	        j++;
	    } else if (tmp[i] == 'I') {
	        ver++;
	        j++;
	        score *= 0.0003;  // 0.0001 + 0.0002;  0.0001: indel rate in normal human, 0.0002: indel error rate in Illumina
	    } else if (tmp[i] == 'D') {
	        ref++;
	        score *= 0.0003; // 0.0001 + 0.0002
	    }
	}
	return score;
}

int matoi(char *str, int start, int end) {
	int i = 0;
	char tmp[SEQ_MAX_LENGTH];
	for (i = 0; i < end - start; i++)
	  tmp[i] = str[start + i];
	tmp[i] = '\0';
	return atoi(tmp);
}

void convertCigarToMatrix(char *cigar, int cigar_size, char * matrix) {
	int i = 0;
	int j = 0;
	int start = 0;
	int size = 0;
	matrix[0] = '\0';
	while (i < cigar_size) {
	    if (cigar[i] >= '0' && cigar[i] <= '9') {
	        start = i;
	        while (cigar[i] >= '0' && cigar[i] <= '9' && i < cigar_size)
	            i++;
	        int value = matoi(cigar, start, i);
	        for (j = 0; j < value; j++) {
	            if (cigar[i] == 'M')
	                matrix[size] = 'M';
	            else if (cigar[i] == 'D')
	                matrix[size] = 'D';
	            else if (cigar[i] == 'I')
	                matrix[size] = 'I';
	            size++;
	        }
	    }
	    i++;
	}
	matrix[size] = '\0';
}

void convertMDToMatrix(char *md, int md_size, char * matrix) {
	  int i = 0;
	  int j = 0;
	  int start = 0;
	  int size = 0;
	  matrix[0] = '\0';
	  while (i < md_size) {
	      if (md[i] >= '0' && md[i] <= '9') {
	          start = i;
	          while (md[i] >= '0' && md[i] <= '9' && i < md_size)
	              i++;
	          int value = matoi(md, start, i);
	          for (j = 0; j < value; j++) {
	              matrix[size] = 'M';
	              size++;
	          }
	          i--;
	      } else if (md[i] == '^') {
	          matrix[size] = 'D';
	          size++;
	      } else {
	          matrix[size] = md[i];
	          size++;
	      }
	      i++;
	  }
	  matrix[size] = '\0';
}

void convertMDCigarToMatrix(char *cigar, int cigar_size, char *md, int md_size,
			    char *matrix) {
	int i = 0;
	int j = 0;
	int size = 0;
	char tmp1[SEQ_MAX_LENGTH];
	char tmp2[SEQ_MAX_LENGTH];
	convertMDToMatrix(md, md_size, tmp2);
	convertCigarToMatrix(cigar, cigar_size, tmp1);
	while (i < strlen(tmp1)) {
	    if (tmp1[i] == 'M') {
	        if (j < strlen(tmp2)) {
	            if (tmp2[j] == 'M') {
	                matrix[size] = 'M';
	                size++;
	            }
	            if (tmp2[j] != 'M') {
	                matrix[size] = tmp2[j];
	                size++;
	            }
	        } else {
	            matrix[size] = 'M';
	            size++;
	        }
	    } else if (tmp1[i] == 'D') {
	        matrix[size] = 'D';
	        size++;
	        j++;
	        matrix[size] = tmp2[j];
	        size++;
	    } else if (tmp1[i] == 'I') {
	        matrix[size] = 'I';
	        size++;
	    }
	    i++;
	    if (j < strlen(tmp2))
	        j++;
	}
	if (strlen(tmp1))
	  matrix[size] = '\0';
}

void convertInsertion(char * in_matrix, char * seq, char *out_matrix) {
	int i = 0;
	int j = 0;
	int size = 0;

	while (i < strlen(in_matrix)) {
	  if (in_matrix[i] == 'M') {
	    out_matrix[size] = 'M';
	    size++;
	    j++;
	  } else if (in_matrix[i] == 'D') {
	    out_matrix[size] = 'D';
	    size++;

	    i++;
	    j++;

	    out_matrix[size] = seq[j];
	    j++;
	    size++;
	  } else if (in_matrix[i] == 'I') {
	    out_matrix[size] = 'I';
	    size++;
	    out_matrix[size] = seq[j];
	    size++;
	    j++;
	  } else {
	    out_matrix[size] = in_matrix[i];
	    size++;
	    j++;
	  }
	  i++;
	}
	out_matrix[size] = '\0';
}

FILE * initPairedEndDiscPP() {
	char fname[FILE_NAME_LENGTH];
	sprintf(fname, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);
	return fileOpen(fname, "w");
}

void finalizePairedEndDiscPP(FILE * _fp_divet) {
	fclose(_fp_divet);
}

void performPairedEndDiscPP(FILE * _fp_divet) {	// out --> _fp_divet
	char tmp_matrix1[SEQ_MAX_LENGTH];
	char tmp_matrix2[SEQ_MAX_LENGTH];

	char matrix1[SEQ_MAX_LENGTH];
	char matrix2[SEQ_MAX_LENGTH];

	char cigar1[MAX_CIGAR_SIZE];
	char editString1[2 * SEQ_MAX_LENGTH];

	char cigar2[MAX_CIGAR_SIZE];
	char editString2[2 * SEQ_MAX_LENGTH];
	char seq1[SEQ_LENGTH + 1];

	char seq2[SEQ_LENGTH + 1];

	char genName[SEQ_LENGTH];
	char fname1[FILE_NAME_LENGTH];
//	char fname2[FILE_NAME_LENGTH];
	char l;
	int l_size;
	int loc1, loc2;
	int err1, err2;
	char dir1, dir2;
	float sc1, sc2, lsc = 0;
	int flag = 0;
	int rNo, lrNo = -1;
	int tmp;
	FILE *in; //, *out;

	sprintf(fname1, "%s__%s__disc", mappingOutputPath, mappingOutput);
//	sprintf(fname2, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);

	in = fileOpen(fname1, "r");
//	out = fileOpen(fname2, "w");

	if (in != NULL) {
	  flag = fread(&rNo, sizeof(int), 1, in);
	} else {
	  flag = 0;
	}

	seq1[SEQ_LENGTH] = '\0';
	seq2[SEQ_LENGTH] = '\0';

	while (flag) {
	  tmp = fread(&l, sizeof(char), 1, in);
	  tmp = fread(genName, sizeof(char), l, in);
	  genName[(int) l] = '\0';
	  tmp = fread(&loc1, sizeof(int), 1, in);
	  tmp = fread(&err1, sizeof(int), 1, in);
	  tmp = fread(&sc1, sizeof(float), 1, in);

	  tmp = fread(&l_size, sizeof(int), 1, in);
		if (l_size > SEQ_LENGTH) 	/* DHL DEBUG */
			break;
	  tmp = fread(cigar1, sizeof(char), l_size, in);
	  cigar1[(int) l_size] = '\0';

	  tmp = fread(&l_size, sizeof(int), 1, in);
	  tmp = fread(editString1, sizeof(char), l_size, in);
		if (l_size > SEQ_LENGTH)	/* DHL DEBUG */
			break;
	  editString1[(int) l_size] = '\0';
	  tmp = fread(&loc2, sizeof(int), 1, in);
	  tmp = fread(&err2, sizeof(int), 1, in);
	  tmp = fread(&sc2, sizeof(float), 1, in);
	  tmp = fread(&l_size, sizeof(int), 1, in);
	  tmp = fread(cigar2, sizeof(char), l_size, in);
		if (l_size > SEQ_LENGTH)	/* DHL DEBUG */
			break;
  	cigar2[(int) l_size] = '\0';
	  tmp = fread(&l_size, sizeof(int), 1, in);
	  tmp = fread(editString2, sizeof(char), l_size, in);
		if (l_size > SEQ_LENGTH)	/* DHL DEBUG */
			break;
  	editString2[(int) l_size] = '\0';

	  convertMDCigarToMatrix(cigar1, strlen(cigar1), editString1,
			   strlen(editString1), tmp_matrix1);
	  convertMDCigarToMatrix(cigar2, strlen(cigar2), editString2,
			   strlen(editString2), tmp_matrix2);
	  /* CHECK FOR SIFAST */
	  /* CALKAN: GO OVER THIS VERY CAREFULLY FOR PE vs MP */
	  if (_msf_readHasConcordantMapping[rNo] == 0 && _msf_discordantMapping[rNo * 2] < maxDiscordantOutput ) {
	      dir1 = dir2 = 'F';
	      strncpy(seq1, _msf_seqList[rNo * 2].seq, SEQ_LENGTH);
	      strncpy(seq2, _msf_seqList[rNo * 2 + 1].seq, SEQ_LENGTH);
	      if (loc1 < 0) {
	          dir1 = 'R';
	          loc1 = -loc1;
	          strncpy(seq1, _msf_seqList[rNo * 2].rseq, SEQ_LENGTH);
	      }
	      if (loc2 < 0) {
	          dir2 = 'R';
	          loc2 = -loc2;
	          strncpy(seq2, _msf_seqList[rNo * 2 + 1].rseq, SEQ_LENGTH);
	      }
	      convertInsertion(tmp_matrix1, seq1, matrix1);
	      convertInsertion(tmp_matrix2, seq2, matrix2);
	      if (rNo != lrNo) {
	          int j;
	          for (j = 0; j < SEQ_LENGTH; j++) {
	              lsc += _msf_seqList[rNo * 2].qual[j]
	                  + _msf_seqList[rNo * 2 + 1].qual[j];
	          }
	          lsc /= 2 * SEQ_LENGTH;
	          lsc -= 33;
	          lrNo = rNo;
	      }
	      char event = '\0';
	      if (dir1 == dir2) {
	          event = 'V';
	      }
	      else {
	          if (pairedEndModePE && loc1 < loc2 && dir1 == 'R' && dir2 == 'F')
	              event = 'E';
	          else if (pairedEndModeMP && loc1 < loc2 && dir1 == 'F' && dir2 == 'R')
	              event = 'E';
	          else if (pairedEndModePE && loc2 < loc1 && dir1 == 'F' && dir2 == 'R')
	              event = 'E';
	          else if (pairedEndModeMP && loc2 < loc1 && dir1 == 'R' && dir2 == 'F')
	              event = 'E';
	          else if (abs(loc2 - loc1) >= maxPairEndedDiscordantDistance)
	              event = 'D';
	          else
	              event = 'I';
	      }
	      _msf_seqList[rNo * 2].hits[0] = 2;
	      fprintf(_fp_divet,
	              "%s\t%s\t%d\t%d\t%c\t=\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%e\n",
	              _msf_seqList[rNo * 2].name, genName, loc1,
	              (loc1 + SEQ_LENGTH - 1), dir1, loc2,
	              (loc2 + SEQ_LENGTH - 1), dir2, event, (err1 + err2),
	              lsc, sc1 * sc2);

	      //	      lsc, sc1 * sc2 * binomial_coefficient(2 * SEQ_LENGTH, (err1 + err2)));
	  }
	  flag = fread(&rNo, sizeof(int), 1, in);
	}

	tmp++;

	fclose(in);
//	fclose(out);

	unlink(fname1);
}

FILE * initOEAReads () {
  FILE * _fp_oea;
	char fname[FILE_NAME_LENGTH];
	sprintf(fname, "%s%s_OEA.sam", mappingOutputPath, mappingOutput);
	_fp_oea = fileOpen(fname, "w");
	return _fp_oea;
}

void finalizeOEAReads(FILE * _fp_oea) {
	fclose(_fp_oea);
}

void performOEAReads(FILE * _fp_oea) {
	FILE* in;
	char genName[SEQ_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	char l = 0;
	int loc1 = 0;
	int err1;
	char d;
	float sc1 = 0;
	int flag = 0;
	int rNo = -1;
	int tmp = 0;
	int cigarSize = 0;
	int mdSize = 0;
	char cigar[SEQ_LENGTH + 1];
	char md[SEQ_LENGTH + 1];
	char *seq1, *seq2, *qual1, *qual2;
	char rqual1[SEQ_LENGTH + 1];

	seq1 = NULL;
	seq2 = NULL;
	qual1 = NULL;
	qual2 = NULL;
	rqual1[0] = '\0';

	SAMheaderTX(_fp_oea, 0);
	in = NULL;
	if (pairedEndDiscordantMode) {
		sprintf(fname2, "%s__%s__oea", mappingOutputPath, mappingOutput);
		in = fileOpen(fname2, "r");
	}

	if (in != NULL) {
		flag = fread(&rNo, sizeof(int), 1, in);
	} else {
		flag = 0;
	}

	while (flag) {
		cigar[0] = '\0';
		md[0] = '\0';

		tmp = fread(&l, sizeof(char), 1, in);
		tmp = fread(genName, sizeof(char), l, in);
		genName[(int) l] = '\0';
		tmp = fread(&loc1, sizeof(int), 1, in);
		tmp = fread(&err1, sizeof(int), 1, in);
		tmp = fread(&sc1, sizeof(float), 1, in);
		tmp = fread(&cigarSize, sizeof(int), 1, in);
		tmp = fread(cigar, sizeof(char), cigarSize, in);
		cigar[cigarSize] = '\0';
		tmp = fread(&mdSize, sizeof(int), 1, in);

		tmp = fread(md, sizeof(char), mdSize, in);
		md[mdSize] = '\0';
		d = 1;

		if (loc1 < 0) {
			d = -1;
			loc1 *= -1;

			seq1 = _msf_seqList[rNo].rseq;
			reverse(_msf_seqList[rNo].qual, rqual1, SEQ_LENGTH);
			rqual1[SEQ_LENGTH] = '\0';
			qual1 = rqual1;
		} else {
			seq1 = _msf_seqList[rNo].seq;
			qual1 = _msf_seqList[rNo].qual;
		}

		if (rNo % 2 == 0) {
			seq2 = _msf_seqList[rNo + 1].seq;
			qual2 = _msf_seqList[rNo + 1].qual;
		} else {
			seq2 = _msf_seqList[rNo - 1].seq;
			qual2 = _msf_seqList[rNo - 1].qual;
		}


		if (_msf_seqHits[rNo] != 0 && _msf_seqHits[rNo] < maxOEAOutput
					&& _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0) {
			_msf_output.POS = loc1;
			_msf_output.MPOS = 0;
			_msf_output.FLAG = (rNo % 2 == 0) ? 1 + 4 + 32 * d + 128 : 1 + 8 + 16 * d + 64;
			_msf_output.ISIZE = 0;
			_msf_output.SEQ = seq1;
			_msf_output.QUAL = qual1;
			_msf_output.QNAME = _msf_seqList[rNo].name;
			_msf_output.RNAME = genName;
			_msf_output.MAPQ = 255;
			_msf_output.CIGAR = cigar;
			_msf_output.MRNAME = "=";
			_msf_output.optSize = 4;
			_msf_output.optFields = _msf_optionalFields;

			_msf_optionalFields[0].tag = "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = err1;

			_msf_optionalFields[1].tag = "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = md;

			//for the OEA reads
			_msf_optionalFields[2].tag = "NS";
			_msf_optionalFields[2].type = 'Z';
			_msf_optionalFields[2].sVal = seq2;

			_msf_optionalFields[3].tag = "NQ";
			_msf_optionalFields[3].type = 'Z';
			_msf_optionalFields[3].sVal = qual2;

			outputSAM(_fp_oea, _msf_output);

			_msf_seqList[rNo].hits[0] = -1;
			_msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
		}
		else if(_msf_seqHits[rNo] != 0 && _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0) {
			_msf_seqList[rNo].hits[0] = -1;
			_msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
		}
		flag = fread(&rNo, sizeof(int), 1, in);
	}
	tmp++;
	fclose(in);
	unlink(fname2);
	//fclose(fp_oea);
}

void outputAllTransChromosomal(int flag) {
	return;
}

void preProcessReadsCG() {

	_msf_sort_seqList = getMem(_msf_seqListSize * sizeof(Pair));

	int i = 0;
	for(i = 0; i < _msf_seqListSize; i++) {
		_msf_sort_seqList[i].hv = hashVal(_msf_seqList[i].seq);
		_msf_sort_seqList[i].readNumber = i;
	}

	qsort(_msf_sort_seqList, _msf_seqListSize, sizeof(Pair), compare);
	_msf_map_sort_seqList = getMem(_msf_seqListSize * sizeof(int));

	for(i = 0; i < _msf_seqListSize; i++) {
		_msf_map_sort_seqList[_msf_seqList[i].readNumber] = i;
	}
}

void initFASTCG(Read *seqList, int seqListSize, int accSeqListSize){
	int i;

	_msf_seqList = seqList;
	_msf_seqListSize = seqListSize;
	_msf_totalSeqListSize = accSeqListSize;

	if (_msf_optionalFields == NULL) {
		_msf_op = getMem(SEQ_LENGTH);

		if (pairedEndMode) {
			_msf_optionalFields = getMem(8*sizeof(OPT_FIELDS));
		}
		else
			_msf_optionalFields = getMem(2*sizeof(OPT_FIELDS));

		for (i = 0; i < 200; i++)
			sprintf(_msf_numbers[i], "%d%c", i, '\0');

		sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
	}


	if (_msf_samplingLocsEnds == NULL) {
		_msf_samplingLocsEnds = getMem(1);

		preProcessReads();

		_msf_oeaMapping = getMem(_msf_seqListSize * sizeof(int));

		for (i = 0; i < _msf_seqListSize; i++)
			_msf_oeaMapping[i] = 0;

		_msf_discordantMapping = getMem(_msf_seqListSize * sizeof(int));

		for (i = 0; i < _msf_seqListSize; i++)
			_msf_discordantMapping[i] = 0;
	}

// Reference Genome Name Update
	if (_msf_refGenName == NULL) {
		_msf_refGenName = getMem(4*SEQ_LENGTH);
	}

	if (_msf_verifiedLocs != NULL) {
		freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength + 1));
		_msf_verifiedLocs = NULL;
	}

	_msf_refGen =  getRefGenome();
	_msf_refGenLength = strlen(_msf_refGen);

	_msf_refGenOffset = getRefGenomeOffset();
	snprintf(_msf_refGenName, 4*SEQ_LENGTH,"%s%c", getRefGenomeName(), '\0');
	_msf_refGenName[strlen(getRefGenomeName())] = '\0';


	_msf_verifiedLocs = (int *) getMem(sizeof(int)*(_msf_refGenLength+1));

	for (i = 0; i <= _msf_refGenLength; i++) {
		_msf_verifiedLocs[i] = _msf_seqListSize*10 + 1;
	}

	if (pairedEndMode && _msf_seqHits == NULL) {
		_msf_mappingInfo  = getMem(seqListSize * sizeof (MappingInfo));

		for (i = 0; i < seqListSize; i++) {
			_msf_mappingInfo[i].next = NULL;
			_msf_mappingInfo[i].size = 0;
		}

		_msf_seqHits = getMem((_msf_seqListSize) * sizeof(int));

		for (i=0; i<_msf_seqListSize; i++) {
			_msf_seqHits[i] = 0;
		}
		_msf_readHasConcordantMapping = getMem(_msf_seqListSize / 2 * sizeof(char));

		for(i = 0; i < _msf_seqListSize/2; i++) {
			_msf_readHasConcordantMapping[i] = 0;
		}
	}

	if (_msf_refGenOffset == 0)
		_msf_refGenBeg = 1;
	else
		_msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2;

	_msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;
}



void mapPairedEndSeqCG() {
	//printf("Before Mapping Memory=%15.2f\n", getMemUsage());
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;

	int key_hash[TEST_KEY_NUM];
	unsigned int *locs = NULL;
	key_struct* key_input = getMem(TEST_KEY_NUM*sizeof(key_struct));

	// First read in forward
	initKeyInput(key_input, key_hash);
	initForwardParam();
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].seq + m*10 + readFormat*5);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
    if (key_input[0].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 1, 2, 3, 0);
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 1, 3, 2, 0);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 2, 3, 1, 0);
    }
    if (key_input[0].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, -1);
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, -1);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, -1);
    }
	}

	initKeyInput(key_input, key_hash);
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].rseq + m*10 + readFormat*5);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
    if (key_input[0].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 1, 2, 3, 0);
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 1, 3, 2, 0);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 2, 3, 1, 0);
    }
    if (key_input[0].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, -1);
      mapPairEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, -1);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, -1);
    }
	}

	// First read in reverse
  initKeyInput(key_input, key_hash);
	initReverseParam();
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].rseq + m*10);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
    if (key_input[0].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, 3);
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, 3);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, 3);
    }
    if (key_input[0].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, -1);
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, -1);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, -1);
    }
	}

  initKeyInput(key_input, key_hash);
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].seq + m*10);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
    if (key_input[0].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, 3);
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, 3);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 1) {
      mapPairEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, 3);
    }
    if (key_input[0].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, -1);
      mapPairEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, -1);
    }
    if (key_input[1].key_entry_size > 0 && readFormat == 0) {
      mapPairEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, -1);
    }
	}
	freeMem(key_input, TEST_KEY_NUM * sizeof(key_struct));

	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	MappingLocations *cur;
	int tmpOut;
	int lmax = 0, rmax = 0;

	sprintf(fname1, "%s__%s__%s__%d__1.tmp", mappingOutputPath, _msf_refGenName,
			mappingOutput, _msf_openFiles);
	sprintf(fname2, "%s__%s__%s__%d__2.tmp", mappingOutputPath, _msf_refGenName,
			mappingOutput, _msf_openFiles);

	FILE* out;
	FILE* out1 = fileOpen(fname1, "w");
	FILE* out2 = fileOpen(fname2, "w");

	_msf_openFiles++;

	for (i = 0; i < _msf_seqListSize; i++) {
		if (i % 2 == 0) {
			out = out1;
			if (lmax < _msf_mappingInfo[i].size) {
				lmax = _msf_mappingInfo[i].size;
			}
		} else {
			out = out2;
			if (rmax < _msf_mappingInfo[i].size) {
				rmax = _msf_mappingInfo[i].size;
			}
		}
		tmpOut = fwrite(&(_msf_mappingInfo[i].size), sizeof(int), 1, out);
		if (_msf_mappingInfo[i].size > 0) {
			cur = _msf_mappingInfo[i].next;
			for (j = 0; j < _msf_mappingInfo[i].size; j++) {
				if (j > 0 && j % MAP_CHUNKS == 0) {
					cur = cur->next;
				}
				if(debugMode && (cur->cigarSize[j % MAP_CHUNKS] > SEQ_LENGTH || cur->mdSize[j % MAP_CHUNKS] > SEQ_LENGTH)) {
					printf("ERROR in %d read size exceeds cigar=%d md =%d cigar=%s md =%s\n", i,  cur->cigarSize[j % MAP_CHUNKS], cur->mdSize[j % MAP_CHUNKS], cur->cigar[j % MAP_CHUNKS], cur->md[j % MAP_CHUNKS]);
				}
				//printf("loc=%d\n", cur->loc[j % MAP_CHUNKS]);
				tmpOut = fwrite(&(cur->loc[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite(&(cur->err[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite(&(cur->cigarSize[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite((cur->cigar[j % MAP_CHUNKS]), sizeof(char), (cur->cigarSize[j % MAP_CHUNKS]), out);
				tmpOut = fwrite(&(cur->mdSize[j % MAP_CHUNKS]), sizeof(int), 1, out);
				tmpOut = fwrite((cur->md[j % MAP_CHUNKS]), sizeof(char), (cur->mdSize[j % MAP_CHUNKS]), out);
			}
		}
	}

	  freeAllMapping();

	_msf_maxLSize += lmax;
	_msf_maxRSize += rmax;
	tmpOut++;


	fclose(out1);
	fclose(out2);

}

// sirfast: searchKeyCG()
int searchKeyCG(int target_coor, unsigned int* entry_coor, int entry_size, int range) {
	if (entry_size <= 0)
		return -1;
	int lower_bound = 1;
	int upper_bound = entry_size;
	int mid = lower_bound + entry_size / 2;

	while (lower_bound < upper_bound) {
		if (entry_coor[mid] == target_coor)
			return entry_coor[mid];
		else if (entry_coor[mid] < target_coor)
			lower_bound = mid + 1;
		else
			upper_bound = mid - 1;
		mid = lower_bound + (upper_bound - lower_bound) / 2;

		if (entry_coor[upper_bound] == target_coor) {
			return entry_coor[upper_bound];
		}
		if (entry_coor[lower_bound] == target_coor) {
			return entry_coor[lower_bound];
		}
		if (entry_coor[mid] == target_coor) {
			return entry_coor[mid];
		}
	}
//	printf ("target_coor: %d, range: %d, entry[mid-1]: %d, entry[mid]: %d, entry[mid+1]:%d \n",
//			target_coor, range, entry_coor[mid-1], entry_coor[mid], entry_coor[mid+1]); // DEBUG
	if (entry_coor[mid] <= (target_coor + range) && entry_coor[mid] >= target_coor) {
		return entry_coor[mid];
	}
	else if (entry_coor[mid+1] <= (target_coor + range) && entry_coor[mid+1] >= target_coor) {
		return entry_coor[mid+1];
	}
	else
		return -1;
}

int verify(int refIndex, int readIndex, char* seq1, int* offset, int variable, int length) {
	int i = 0;
	int j = 0;
	int errMin = 10;
	//fprintf(stdout, "variable: %d\n", variable);

	for (i = 0; i < variable; i++) {
		int errCnt = 0;
		char* ref = _msf_refGen + refIndex + i - 1;
		char* seq = seq1 + readIndex;
		for (j = 0; j < length; j++) {
			//fprintf(stdout, "[%c:%c]\n", *ref, *seq);
			errCnt += *ref != *seq;
			seq++;
			ref++;
		}
		//for (j = 0; j < length; j++)
		//	fprintf(stdout, "%c", *(_msf_refGen + refIndex + i + j));
		//fprintf(stdout,"\n");
		//for (j = 0; j < length; j++)
		//	fprintf(stdout, "%c", *(seq1 + readIndex + j));
		//fprintf(stdout,"\n");
		//fprintf(stdout, "  ERROR : %d", errCnt);
		//fprintf(stdout, "  OFFSET: %d\n", i);
		if (errCnt < errMin) {
			errMin = errCnt;
			offset[0] = i;
		}
	}
	return errMin;
}

void initKeyInput(key_struct* key_input, int* key_hash) {
	int i = 0;
	for(i = 0; i < TEST_KEY_NUM; i++) {
		key_input[i].key_entry = NULL;
		key_input[i].key_locs = NULL;
		key_input[i].key_number = 0;
		key_input[i].key_entry_size = 0;
		key_hash[i] = 0;
	}
}

void initForwardParam() {
	if (readFormat == 1) {
		initial_gap[0] = -1;
		initial_gap[1] =  0;
		initial_gap[2] =  5;
		expected_gap[0] = 2;
		expected_gap[1] = 2;
		expected_gap[2] = 2;
		key_length[0] = 5;
		key_length[1] = 10;
		key_length[2] = 10;
		key_length[3] = 10;
	}
	else {
		initial_gap[0] = -7;
		initial_gap[1] = -1;
		expected_gap[0] = 6;
		expected_gap[1] = 0;
		//expected_gap[2] = 0;
		key_length[0] = 10;
		key_length[1] = 10;
		key_length[2] = 10;
	}
}

void initReverseParam() {
	if (readFormat == 1) {
		initial_gap[0] =  5;
		initial_gap[1] =  0;
		initial_gap[2] = -3;
		expected_gap[0] = 2;
		expected_gap[1] = 2;
		expected_gap[2] = 2;
		key_length[0] = 10;
		key_length[1] = 10;
		key_length[2] = 10;
		key_length[3] = 5;
	}
	else {
		initial_gap[0] = -1;
		initial_gap[1] = -7;
		expected_gap[0] = 0;
		expected_gap[1] = 6;
		//expected_gap[2] = 6;
		key_length[0] = 10;
		key_length[1] = 10;
		key_length[2] = 10;
	}
}

// sirfast: mapAllSingleEndSeqCG()
void mapAllSingleEndSeqCG() {
	int i = 0;
	int k = 0;
	int m = 0;
	int key_hash[TEST_KEY_NUM];
	unsigned int *locs = NULL;
	key_struct* key_input = getMem(TEST_KEY_NUM*sizeof(key_struct));
	//readFormat = 0;

	// Forward
	initKeyInput(key_input, key_hash);
	initForwardParam();
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].seq + m*10 + readFormat*5);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 1, 2, 3, 0);
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 1, 3, 2, 0);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 2, 3, 1, 0);
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, -1);
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, -1);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, -1);
		}
	}

	initKeyInput(key_input, key_hash);
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].rseq + m*10 + readFormat*5);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number	= m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 1, 2, 3, 0);
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 1, 3, 2, 0);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 2, 3, 1, 0);
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, -1);
			mapSingleEndSeqCG_forward(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, -1);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_forward(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, -1);
		}
	}

	//Reverse
	initKeyInput(key_input, key_hash);
	initReverseParam();

	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].rseq + m*10);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number	= m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, 3);
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, 3);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, 3);
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 1, 2, -1);
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 1, 0, 2, 1, -1);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 1, 1, 2, 0, -1);
		}
	}

	initKeyInput(key_input, key_hash);
	for(i = 0; i < _msf_seqListSize; i++) {
		k = _msf_sort_seqList[i].readNumber;
		for (m = 0; m < TEST_KEY_NUM; m++) {
			key_hash[m] = hashVal(_msf_seqList[k].seq + m*10);
			locs = getCandidates(key_hash[m]);
			key_input[m].key_number = m;
			key_input[m].key_entry = locs;
			key_input[m].key_locs = locs + 1;
			key_input[m].key_entry_size = -1;
			if (locs != NULL)
				key_input[m].key_entry_size = locs[0];
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, 3);
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, 3);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 1) {
			mapSingleEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, 3);
		}
		if (key_input[0].key_entry_size > 0 && readFormat == 0) {
	  	mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 1, 2, -1);
			mapSingleEndSeqCG_reverse(key_input[0].key_locs, key_input[0].key_entry_size, k, key_input[0].key_number, 0, key_input, 0, 0, 2, 1, -1);
		}
		if (key_input[1].key_entry_size > 0 && readFormat == 0) {
			mapSingleEndSeqCG_reverse(key_input[1].key_locs, key_input[1].key_entry_size, k, key_input[1].key_number, 1, key_input, 0, 1, 2, 0, -1);
		}
	}

	freeMem(key_input, TEST_KEY_NUM * sizeof(key_struct));
	return ;
}

void printoutCG_forward(int genLoc, int* af_offset, char* seq1, int error) {
	char* ref;
	int n = 0;
	int m = 0;
	int ix = 0;
	int iy = 0;
	int is = 0;
//	int space[4];
	int length[4];

	ref = _msf_refGen + genLoc;
	for (ix = 0; ix < 4; ix++) {
//		space[ix] = EXPECTED_GAP + 1;
		length[ix] = KEY_LENGTH;
	}
	length[0] = 5 + af_offset[0];

	n = 0;
	fprintf(stdout, "\nForward Mapping");
	fprintf(stdout, "\nReference: ");

	for (n = 0; n < 4; n++) {
		for (ix = is; ix < is + length[n]; ix++) {
			fprintf(stdout, "%c", *(ref + ix));
		}
		is = ix;
		fprintf(stdout, " ");
		for (m = 0; m < -af_offset[n]; m++) {
			fprintf(stdout, " ");
		}
		if (n < 3) {
			for (ix = is; ix < is + af_offset[n]; ix++) {
				fprintf(stdout, "%c", *(ref + ix));
			}
		}
		is = ix;
		fprintf(stdout, " ");
	}
	fprintf(stdout,"\n");

	n = 0;
	fprintf(stdout, "Fragment : ");
	for(ix = 0; ix < 35; ix++) {
		fprintf(stdout, "%c", seq1[ix]);
		if (ix == 4 ) {
			for (iy = 0; iy < 2; iy++) {
				fprintf(stdout," ");
			}
			n = n + 1;
		}
		if (ix == 14 || ix == 24) {
			for (iy = 0; iy < 2 + af_offset[n]; iy++) {
				fprintf(stdout," ");
			}
			n = n + 1;
		}
	}
	fprintf(stdout,"\n");
	fprintf(stdout, "---> OFFSET:%d:%d:%d \n", af_offset[0], af_offset[1], af_offset[2]);
	fprintf(stdout, "---> LOC=%d refGENEOFFSET=%d\n", genLoc , _msf_refGenOffset);
	fprintf(stdout, "---> MAPLOC:%d \n", genLoc + _msf_refGenOffset);
	fprintf(stdout, "---> ERROR :%d \n\n", error);
	fprintf(stdout, "%s\t%s\t%d\t%d5M%dB10M%dN10M%dN10\n", seq1, _msf_refGenName, genLoc + _msf_refGenOffset, error , -af_offset[0], af_offset[1], af_offset[2]);
}

void printoutCG_reverse(int genLoc, int * af_offset, char * seq1, int error) {
	char * ref;
	int n = 0;
	int m = 0;
	int ix = 0;
	int iy = 0;
	int is = 0;
//	int space[4];
	int length[4];

	ref = _msf_refGen + genLoc;
	for (ix = 0; ix < 4; ix++) {
//		space[ix] = EXPECTED_GAP + 1;
		length[ix] = KEY_LENGTH;
	}
	length[3] = 5 + af_offset[2];

	n = 0;
	fprintf(stdout, "\nReverse Mapping");
	fprintf(stdout, "\nReference: ");
	for (n = 0; n < 4; n++) {
		for (ix = is; ix < is + length[n]; ix++) {
			fprintf(stdout, "%c", *(ref + ix));
		}
		is = ix;
		fprintf(stdout, " ");
		for (m = 0; m < -af_offset[n]; m++) {
			fprintf(stdout, " ");
		}
		if (n < 3) {
			for (ix = is; ix < is + af_offset[n]; ix++) {
				fprintf(stdout, "%c", *(ref + ix));
			}
		}
		is = ix;
		fprintf(stdout, " ");
	}
	fprintf(stdout,"\n");

	n = 0;
	fprintf(stdout, "Fragment : ");
	for(ix = 0; ix < 35; ix++) {
		fprintf(stdout, "%c", seq1[ix]);
		if (ix == 29 ) {
			for (iy = 0; iy < 2; iy++) {
				fprintf(stdout," ");
			}
			n = n + 1;
		}
		if (ix == 9 || ix == 19) {
			for (iy = 0; iy < 2 + af_offset[n]; iy++) {
				fprintf(stdout," ");
			}
			n = n + 1;
		}
	}
	fprintf(stdout,"\n");
	fprintf(stdout, "---> OFFSET:%d:%d:%d \n", af_offset[0], af_offset[1], af_offset[2]);
	fprintf(stdout, "---> LOC=%d refGENEOFFSET=%d\n", genLoc , _msf_refGenOffset);
	fprintf(stdout, "---> MAPLOC:%d \n", genLoc + _msf_refGenOffset);
	fprintf(stdout, "---> ERROR :%d \n\n", error);
	fprintf(stdout, "%s\t%s\t%d\t%d5M%dB10M%dN10M%dN10\n", seq1, _msf_refGenName, genLoc + _msf_refGenOffset, error , -af_offset[0], af_offset[1], af_offset[2]);
}

void mapPairEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index,
					key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext) {

	char matrix[200];
	char editString[200];
	char cigar[MAX_CIGAR_SIZE];

	int r = readNumber;
	int d = (direction==1?-1:1);

	int readId = 2*readNumber + direction;

	char *_tmpSeq;
	char rqual[SEQ_LENGTH+1];
	rqual[SEQ_LENGTH]='\0';

	if (direction) {
		reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
		_tmpSeq = _msf_seqList[readNumber].rseq;
	}
	else {
		_tmpSeq = _msf_seqList[readNumber].seq;
	}

	int i = 0;
	int j = 0;
	int genLoc = 0;
	int *locs = (int *) l1;

	for (j = 0; j < s1; j++) {
		genLoc = locs[j];
		int af_offset[4];
		int err = -1;

		if (readFormat == 1) {
      if (key_input[index].key_number == 0) {
				if (genLoc - 1			< _msf_refGenBeg
				 || genLoc - 1 + 30 + 8 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10	- 7	< _msf_refGenBeg
				 || genLoc - 1 + 25 - 1 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 10 - 6] == readId ) {
					continue;
				}
			}
		}
		else {
      if (key_input[index].key_number == 0) {
				if (genLoc - 1 < _msf_refGenBeg
				 || genLoc - 1 + 30 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10	< _msf_refGenBeg
				 || genLoc - 1 + 20 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 10 - 1] == readId ) {
					continue;
				}
			}
		}

	  err = verifySingleEndSeqCG_backward(&genLoc, af_offset, _tmpSeq, key_input, index, num_base, num_key, num_ver, num_ext);

		if (err <= errThreshold && err >= 0 && _msf_verifiedLocs[genLoc] != readId) {
			generateAlignmentMatrxCG_backward(genLoc, af_offset, _tmpSeq, err, matrix);
			generateSNPSAM(matrix, strlen(matrix), editString);
			if (readFormat == 1)
				sprintf(cigar, "5M%dS10M%dN10M%dN10M", -af_offset[2], af_offset[1], af_offset[0]);
			else
				sprintf(cigar, "10M%dS10M%dS10M", -af_offset[1], -af_offset[0]);

		}
		else {
			err = -1;
		}


		//##### mrfast code #####

		if(err != -1 && !bestMode) {
			int offset_range = 3;
			for(i = -offset_range ; i <= offset_range ; i++) {
				if(genLoc + i >= _msf_refGenBeg && genLoc + i <= _msf_refGenEnd) {
					_msf_verifiedLocs[genLoc + i] = readId;
				}
			}

			/* calkan counter */
			mappingCnt++;
			MappingLocations *parent = NULL;
			MappingLocations *child = _msf_mappingInfo[r].next;

			for (i = 0; i < (_msf_mappingInfo[r].size / MAP_CHUNKS); i++) {
				parent = child;
				child = child->next;
			}

			if (child == NULL) {
				MappingLocations *tmp = getMem(sizeof(MappingLocations));
				tmp->next = NULL;
				tmp->loc[0] = (genLoc+_msf_refGenOffset) * d; // d is required
				tmp->err[0] = err;
				tmp->cigarSize[0] = strlen(cigar);
				sprintf(tmp->cigar[0], "%s", cigar);
				tmp->mdSize[0] = strlen(editString);	// DHL
				sprintf(tmp->md[0], "%s", editString);

				if (parent == NULL)
					_msf_mappingInfo[r].next = tmp;
				else
					parent->next = tmp;
			} else {
				if (strlen(cigar) > SEQ_LENGTH
						|| strlen(editString) > SEQ_LENGTH) {
					printf(
							"ERROR in %d read size(After mapping) exceeds cigar=%d md =%d cigar=%s md =%s\n",
							r, (int) strlen(cigar), (int) strlen(editString),
							cigar, editString);
				}
				child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = (genLoc+_msf_refGenOffset) * d;
				child->err[_msf_mappingInfo[r].size % MAP_CHUNKS] = err;
				child->cigarSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(cigar);
				sprintf(child->cigar[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", cigar);
				child->mdSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(editString);
				sprintf(child->md[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", editString);
			}
			_msf_mappingInfo[r].size++;

		}
	}
}

void freeAllMapping()
{
	  int i = 0;
	  int j = 0;

	  MappingLocations *prev;
	  MappingLocations *cur;

	  //printf("Before Cleaning Memory=%15.2f\n", getMemUsage());

	  for(i = 0; i < _msf_seqListSize; i++) {
	       if (_msf_mappingInfo[i].size > 0) {
	          cur = _msf_mappingInfo[i].next;
	          for(j = 0; j < _msf_mappingInfo[i].size; j++) {
	              if(j>0 && j % MAP_CHUNKS == 0) {
	                  prev = cur;
	                  cur = cur->next;
	                  freeMem(prev, sizeof(MappingLocations));
	              }
	          }
	          if(cur != NULL)
	              freeMem(cur, sizeof(MappingLocations));
	       }
	       _msf_mappingInfo[i].next = NULL;
	       _msf_mappingInfo[i].size = 0;
	  }
	 // printf("After Cleaning Memory=%15.2f\n", getMemUsage());
	 // printf("Cleaning All\n");
}

void mapPairEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index,
				key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext) {

	char matrix[MAX_CIGAR_SIZE];
	char editString[MAX_CIGAR_SIZE];
	char cigar[MAX_CIGAR_SIZE];

	int r = readNumber;
	int d = (direction==1?-1:1);

	int readId = 2*readNumber+direction;
	char *_tmpSeq;
	char rqual[SEQ_LENGTH+1];
	rqual[SEQ_LENGTH]='\0';

	int i = 0;

	if (direction) {
		reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
		_tmpSeq = _msf_seqList[readNumber].rseq;
	}
	else {
		_tmpSeq = _msf_seqList[readNumber].seq;
	}

	int j = 0;
	int genLoc = 0;
	int *locs = (int *) l1;

	for (j = 0; j < s1; j++) {
		genLoc = locs[j];
		int af_offset[4];
		int err = -1;

		if (readFormat == 1) {
			if (key_input[index].key_number == 0) {
				if (genLoc - 1 - 5 + 1  < _msf_refGenBeg
				 || genLoc - 1 + 30 + 9 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 5] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 15 - 1 < _msf_refGenBeg
				 || genLoc - 1 + 20 + 7 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 15] == readId ) {
					continue;
				}
			}
		}
		else {
			if (key_input[index].key_number == 0) {
				if (genLoc - 1 + 1  < _msf_refGenBeg
				 || genLoc - 1 + 30 + 9 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10 - 1 < _msf_refGenBeg
				 || genLoc - 1 + 20 + 7 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 1] == readId ) {
					continue;
				}
			}
		}

		err = verifySingleEndSeqCG_forward(&genLoc, af_offset, _tmpSeq, key_input, index, num_base, num_key, num_ver, num_ext);

		if (err <= errThreshold && err >= 0 && _msf_verifiedLocs[genLoc] != readId) {
			generateAlignmentMatrxCG_forward(genLoc, af_offset, _tmpSeq, err, matrix);
			generateSNPSAM(matrix, strlen(matrix), editString);
			if (readFormat == 1)
				sprintf(cigar, "5M%dS10M%dN10M%dN10M", -af_offset[0], af_offset[1], af_offset[2]);
			else
				sprintf(cigar, "10M%dS10M%dS10M", -af_offset[0], -af_offset[1]);
		}
		else {
			err = -1;
		}

		//##### mrfast code #####

		if(err != -1 && !bestMode) {
			int offset_range = 3;
			for(i = -offset_range ; i <= offset_range ; i++) {
				if(genLoc + i >= _msf_refGenBeg && genLoc + i <= _msf_refGenEnd) {
					_msf_verifiedLocs[genLoc + i] = readId;
				}
			}

			/* calkan counter */
			mappingCnt++;
			MappingLocations *parent = NULL;
			MappingLocations *child = _msf_mappingInfo[r].next;

			for (i = 0; i < (_msf_mappingInfo[r].size / MAP_CHUNKS); i++) {
				parent = child;
				child = child->next;
			}

			if (child == NULL) {
				//printf("NULL child\n");
				MappingLocations *tmp = getMem(sizeof(MappingLocations));
				tmp->next = NULL;
				tmp->loc[0] = (genLoc+_msf_refGenOffset) * d; // d is required
				tmp->err[0] = err;
				tmp->cigarSize[0] = strlen(cigar);
				sprintf(tmp->cigar[0], "%s", cigar);
				tmp->mdSize[0] = strlen(editString);
				sprintf(tmp->md[0], "%s", editString);
				if (parent == NULL)
					_msf_mappingInfo[r].next = tmp;
				else
					parent->next = tmp;
			} else {
				//printf("NOT NULL CHILD\n");
				if (strlen(cigar) > SEQ_LENGTH || strlen(editString) > SEQ_LENGTH) {
					printf(
							"ERROR in %d read size(After mapping) exceeds cigar=%d md =%d cigar=%s md =%s\n",
							r, (int) strlen(cigar), (int) strlen(editString),
							cigar, editString);
				}
				child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = (genLoc + _msf_refGenOffset) * d;
				child->err[_msf_mappingInfo[r].size % MAP_CHUNKS] = err;
				child->cigarSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(cigar);
				sprintf(child->cigar[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", cigar);
				child->mdSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(editString);
				sprintf(child->md[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", editString);
			}
			_msf_mappingInfo[r].size++;
		}
	}

}

void generateAlignmentMatrxCG_backward(int genLoc, int * af_offset, char * seq, int error, char * matrix) {
	char * ref;
	int ix = 0;
	ref = _msf_refGen + genLoc;

	if (readFormat == 1) {	// Format: 5-10-10-10
		for(ix = 0; ix < 10; ix++) {
			if(ref[ix] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix];
		}
		for(ix = 10; ix < 20; ix++) {
			if(ref[ix+af_offset[0]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]];
		}
		for(ix = 20; ix < 30; ix++) {
			if(ref[ix + af_offset[0] + af_offset[1]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix + af_offset[0] + af_offset[1]];
		}
		for(ix = 30; ix < 35; ix++) {
			if(ref[ix + af_offset[0] + af_offset[1] + af_offset[2]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix + af_offset[0] + af_offset[1] + af_offset[2]];
		}
	}
	else {	// Format: 10-9-N-10
		for(ix = 0; ix < 10; ix++) {
			if(ref[ix] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix];
		}
		for(ix = 10; ix < 20; ix++) {
			if(ref[ix+af_offset[0]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]];
		}
		for(ix = 20; ix < 30; ix++) {
			if(ref[ix+af_offset[0]+af_offset[1]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]+af_offset[1]];
		}
	}
	matrix[ix] = '\0';
}

void generateAlignmentMatrxCG_forward(int genLoc, int * af_offset, char * seq, int error, char * matrix) {
	char * ref;
	int ix = 0;
	ref = _msf_refGen + genLoc;

	if (readFormat == 1) {	// Format: 5-10-10-10
		for(ix = 0; ix < 5; ix++) {
			if(ref[ix] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix];
		}
		for(ix = 5; ix < 15; ix++) {
			if(ref[ix+af_offset[0]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]];
		}
		for(ix = 15; ix < 25; ix++) {
			if(ref[ix + af_offset[0] + af_offset[1]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]+af_offset[1]];
		}
		for(ix = 25; ix < 35; ix++) {
			if(ref[ix + af_offset[0] + af_offset[1] + af_offset[2]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix + af_offset[0] + af_offset[1] + af_offset[2]];
		}
	}
	else {	// Format: 10-9-N-10
		for(ix = 0; ix < 10; ix++) {
			if(ref[ix] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix];
		}
		for(ix = 10; ix < 20; ix++) {
			if(ref[ix+af_offset[0]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix+af_offset[0]];
		}
		for(ix = 20; ix < 30; ix++) {
			if(ref[ix + af_offset[0] + af_offset[1]] == seq[ix])
				matrix[ix] = 'M';
			else
				matrix[ix] = ref[ix + af_offset[0] + af_offset[1]];
		}
	}
	matrix[ix] = '\0';
}

int verifySingleEndSeqCG_forward(int* locs, int* af_offset, char* seq1, key_struct* key_input,
																 int index, int num_base, int num_key, int num_ver, int num_ext) {
	int i = 0;
	int err = -1;
	int offset[1] = {0};
	int genLoc = *locs;

	// Status 1: Initial Variable Setup
	af_offset[0] = initial_gap[0];
	af_offset[1] = initial_gap[1];
	af_offset[2] = initial_gap[2];

	// Forward

	// ========================================================== //
	//                  BASE            KEY            VERIFY
	// (seq0) -------- (seq1) -------- (seq2) -------- (seq3)
	//  len0            len1            len2            len3
	//          gap0            gap1            gap2
	// ========================================================== //
	//                 VERIFY           BASE            KEY
	// (seq0) -------- (seq1) -------- (seq2) -------- (seq3)
	//  len0            len1            len2            len3
	//          gap0            gap1            gap2
	// ========================================================== //


	// Status 2: Search Key
	// -- expected key location & expected gap variable
	int gap_key = 0;
	int exp_key = genLoc;
	for (i = num_base; i < num_key; i++) {
		exp_key += key_length[i] + initial_gap[i];
		gap_key += expected_gap[i];
	}

	int loc_key = 0;
	if (readFormat == 1)
		loc_key = searchKeyCG(exp_key, key_input[num_key-1].key_entry, key_input[num_key-1].key_entry_size, gap_key);
	else
		loc_key = searchKeyCG(exp_key, key_input[num_key].key_entry, key_input[num_key].key_entry_size, gap_key);

	if (loc_key == -1)
		return -1;

	af_offset[num_key - 1] = loc_key - exp_key + af_offset[num_key - 1];

	// Status 3: Verify
	// -- expected verification location
	int seq_var = 0;
	for (i = 0; i < num_ver; i++) {
		seq_var += key_length[i];
	}

	int exp_var = genLoc;
	for (i = num_base; i < num_ver; i++) {
		exp_var += af_offset[i] + key_length[i];
	}
	for (i = num_ver; i < num_base; i++) {
		exp_var -= af_offset[i] + expected_gap[i] + key_length[i];
	}

	if (num_key - num_base == 2 && readFormat == 1) {
		int tmp_offset = af_offset[num_key - 1];
		int threshold = 0;
		for (i = num_base; i < num_key; i++) {
			threshold += initial_gap[i] + expected_gap[i]/2;
		}
		int tmp_variable = (threshold > tmp_offset ? threshold - tmp_offset : tmp_offset - threshold);

		if (tmp_offset <= threshold) {
			err = verify(exp_var, seq_var, seq1, offset, expected_gap[num_ver] + 1 - tmp_variable, key_length[num_ver]);
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + *offset;
			af_offset[num_ver] = tmp_offset - af_offset[num_ver - 1];
		}
		else {
			err = verify(exp_var + tmp_variable, seq_var, seq1, offset, expected_gap[num_ver] + 1 - tmp_variable, key_length[num_ver]);
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + tmp_variable + *offset;
			af_offset[num_ver] = tmp_offset - af_offset[num_ver - 1];
		}
	}
	else if (num_key - num_base == 2 && readFormat != 1) {
		af_offset[0] = af_offset[0] + af_offset[1] - initial_gap[1];
		af_offset[1] = initial_gap[1];
		exp_var = genLoc + af_offset[0] + key_length[0];
		err = verify(exp_var, seq_var, seq1, offset, 1, key_length[num_ver]);
	}
	else {
		err = verify(exp_var, seq_var, seq1, offset, expected_gap[num_ver] + 1, key_length[num_ver]);
		if (num_ver > num_base) {
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + *offset;
		}
		else {
			af_offset[num_ver] = initial_gap[num_ver] + expected_gap[num_ver] - *offset;
		}
	}

	if (readFormat == 1) {
		int exp_ext = genLoc;
		for (i = 0; i < num_base; i++) {
			exp_ext -= key_length[i] + af_offset[i];
		}
		int seq_ext = 0;
		for (i = 0; i < num_ext; i++) {
			seq_ext += key_length[i];
		}
		err += verify(exp_ext, seq_ext, seq1, offset, expected_gap[num_ext] + 1, key_length[num_ext]);
		af_offset[0] = initial_gap[0] - *offset;
	}

	if (err > errThreshold)
		return -1;

	for (i = 0; i < num_base; i++)
		genLoc -= key_length[i] + af_offset[i];
	genLoc -= 1;
	*locs = genLoc;

	//printoutCG_forward(genLoc, af_offset, seq1, err);
	return err;

}

int verifySingleEndSeqCG_backward(int* locs, int* af_offset, char* seq1, key_struct* key_input,
																	int index, int num_base, int num_key, int num_ver, int num_ext) {
	//printf("num_base: %d, num_key: %d, num_ver: %d, num_ext: %d\n", num_base, num_key, num_ver, num_ext);
	int i = 0;
	int err = -1;
	int offset[1] = {0};
	int genLoc = *locs;

	// Status 1: Initial Variable Setup
	af_offset[0] = initial_gap[0];
	af_offset[1] = initial_gap[1];
	af_offset[2] = initial_gap[2];
	//printf("initial_gap0: %d, initial_gap1: %d, initial_gap: %d\n", initial_gap[0], initial_gap[1], initial_gap[2]);

	// Reverse
	// ========================================================== //
	//  BASE            KEY             VERIFY
	// (seq0) -------- (seq1) -------- (seq2) -------- (seq3)
	//  len0            len1            len2            len3
	//          gap0            gap1            gap2
	// ========================================================== //

	// Status 2: Search Key
	// -- expected key location & expected gap variable
	int gap_key = 0;
	int exp_key = genLoc;
	for (i = num_base; i < num_key; i++) {
		exp_key += key_length[i] + initial_gap[i];
		gap_key += expected_gap[i];
	}

	int loc_key = searchKeyCG(exp_key, key_input[num_key].key_entry, key_input[num_key].key_entry_size, gap_key);
	//printf("base: %d, exp_key: %d, loc_key: %d gap_key: %d\n", genLoc, exp_key, loc_key, gap_key);
	if (loc_key == -1)
		return -1;

	af_offset[num_key-1] = loc_key - exp_key + af_offset[num_key-1];

	// Status 3: Verify
	// -- expected verification location
	int seq_var = 0;
	for (i = 0; i < num_ver; i++) {
		seq_var += key_length[i];
	}

	int exp_var = genLoc;
	for (i = num_base; i < num_ver; i++) {
		exp_var += af_offset[i] + key_length[i];
	}
	for (i = num_ver; i < num_base; i++) {
		exp_var -= af_offset[i] + expected_gap[i] + key_length[i];
	}

	if (num_key - num_base == 2 && readFormat == 1) {
		int tmp_offset = af_offset[num_key - 1];
		int threshold = 0;
		for (i = num_base; i < num_key; i++) {
			threshold += initial_gap[i] + expected_gap[i]/2;
		}
		int tmp_variable = (threshold > tmp_offset ? threshold - tmp_offset : tmp_offset - threshold);

		if (tmp_offset <= threshold) {
			err = verify(exp_var, seq_var, seq1, offset, expected_gap[num_ver] + 1 - tmp_variable, key_length[num_ver]);
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + *offset;
			af_offset[num_ver] = tmp_offset - af_offset[num_ver - 1];
		}
		else {
			err = verify(exp_var + tmp_variable, seq_var, seq1, offset, expected_gap[num_ver] + 1 - tmp_variable, key_length[num_ver]);
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + tmp_variable + *offset;
			af_offset[num_ver] = tmp_offset - af_offset[num_ver - 1];
		}
	}
	else if (num_key - num_base == 2 && readFormat != 1) {
		exp_var = genLoc + key_length[0] + initial_gap[0];
		err = verify(exp_var, seq_var, seq1, offset, 1, key_length[num_ver]);
	}
	else {
		if (num_ver > num_base) {
			err = verify(exp_var, seq_var, seq1, offset, expected_gap[num_ver-1] + 1, key_length[num_ver]);
			af_offset[num_ver - 1] = initial_gap[num_ver - 1] + *offset;
		}
		else {
			err = verify(exp_var, seq_var, seq1, offset, expected_gap[num_ver] + 1, key_length[num_ver]);
			af_offset[num_ver] = initial_gap[num_ver] + expected_gap[num_ver] - *offset;
		}
	}

	if (readFormat == 1) {
		int exp_ext = genLoc;
		for (i = num_base; i < num_ext; i++)
			exp_ext += key_length[i] + af_offset[i];
		int seq_ext = 0;
		for (i = 0; i < num_ext; i++) {
			seq_ext += key_length[i];
		}
		err += verify(exp_ext, seq_ext, seq1, offset, expected_gap[num_ext-1] + 1, key_length[num_ext]);
		af_offset[num_ext] = initial_gap[num_ext] - *offset;
	}

	if (err > errThreshold)
		return -1;

	for (i = 0; i < num_base; i++)
		genLoc -= key_length[i] + af_offset[i];
	genLoc -= 1;

	*locs = genLoc;

	//printoutCG_reverse(genLoc, af_offset, seq1, err);
	return err;
}

// sirfast: mapSingleEndSeqCG_forward
// first_mate 0 or 1, 0 the read is the first part and 1 is the second part.
void mapSingleEndSeqCG_forward(unsigned int *l1, int s1, int readNumber, int readSegment, int index,
						key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext) {

	int i = 0;
	int j = 0;
	char matrix[MAX_CIGAR_SIZE];
	char editString[MAX_CIGAR_SIZE];
	char cigar[MAX_CIGAR_SIZE];
	int readId = 2*readNumber+direction;
	char *_tmpSeq, *_tmpQual;
	char rqual[SEQ_LENGTH+1];
	rqual[SEQ_LENGTH]='\0';
	int genLoc = 0;
	int *locs = (int *) l1;

	if (direction) {
		reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
		_tmpQual = rqual;
		_tmpSeq = _msf_seqList[readNumber].rseq;
	}
	else {
		_tmpQual = _msf_seqList[readNumber].qual;
		_tmpSeq = _msf_seqList[readNumber].seq;
	}

	for (j = 0; j < s1; j++) {
		genLoc = locs[j];
		int af_offset[4];
		int err = -1;
		if (readFormat == 1) {
			if (key_input[index].key_number == 0) {
				if (genLoc - 1 - 5 + 1  < _msf_refGenBeg
				 || genLoc - 1 + 30 + 9 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 5] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 15 - 1 < _msf_refGenBeg
				 || genLoc - 1 + 20 + 7 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 15] == readId ) {
					continue;
				}
			}
		}
		else {
			if (key_input[index].key_number == 0) {
				if (genLoc - 1 + 1  < _msf_refGenBeg
				 || genLoc - 1 + 30 + 9 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10 - 1 < _msf_refGenBeg
				 || genLoc - 1 + 20 + 7 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 1] == readId ) {
					continue;
				}
			}
		}

		err = verifySingleEndSeqCG_forward(&genLoc, af_offset, _tmpSeq, key_input, index, num_base, num_key, num_ver, num_ext);

		if (err <= errThreshold && err >= 0 && _msf_verifiedLocs[genLoc] != readId ) {
			generateAlignmentMatrxCG_forward(genLoc, af_offset, _tmpSeq, err, matrix);
			generateSNPSAM(matrix, strlen(matrix), editString);
			if (readFormat == 1)
	      sprintf(cigar, "5M%dS10M%dN10M%dN10M", -af_offset[0], af_offset[1], af_offset[2]);
			else
	      sprintf(cigar, "10M%dS10M%dS10M", -af_offset[0], -af_offset[1]);
	  }
		else {
			err = -1;
		}

		if(err != -1 && !bestMode) {
			mappingCnt++;
			int offset_range = 3;
			for(i = -offset_range ; i <= offset_range ; i++) {
				if(genLoc + i >= _msf_refGenBeg && genLoc + i <= _msf_refGenEnd) {
					_msf_verifiedLocs[genLoc + i] = readId;
				}
			}

			_msf_seqList[readNumber].hits[0]++;
			_msf_output.QNAME		= _msf_seqList[readNumber].name;
			_msf_output.FLAG		= 16 * direction;
			_msf_output.RNAME		= _msf_refGenName;
			_msf_output.POS			= genLoc + _msf_refGenOffset;
			_msf_output.MAPQ		= 255;
			_msf_output.CIGAR		= cigar;
			_msf_output.MRNAME		= "*";
			_msf_output.MPOS		= 0;
			_msf_output.ISIZE		= 0;
			_msf_output.SEQ			= _tmpSeq;
			_msf_output.QUAL		= _tmpQual;
			_msf_output.optSize		= 2;
			_msf_output.optFields	= _msf_optionalFields;
			_msf_optionalFields[0].tag = "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = err;
			_msf_optionalFields[1].tag = "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = editString;

			output(_msf_output);

			if (_msf_seqList[readNumber].hits[0] == 1) {
				mappedSeqCnt++;
			}
			if ( maxHits == 0 ) {
				_msf_seqList[readNumber].hits[0] = 2;
			}

			if ( maxHits!=0 && _msf_seqList[readNumber].hits[0] == maxHits) {
				completedSeqCnt++;
				break;
			}
		}
	}
}

// sirfast: mapSingleEndSeqCG_reverse
// first_mate 0 or 1, 0 the read is the first part and 1 is the second part.
void mapSingleEndSeqCG_reverse(unsigned int *l1, int s1, int readNumber, int readSegment, int index,
						key_struct* key_input, int direction, int num_base, int num_key, int num_ver, int num_ext) {
	char matrix[200];
	char editString[200];
	char cigar[MAX_CIGAR_SIZE];

	int readId = 2 * readNumber + direction;

	char *_tmpSeq, *_tmpQual;
	char rqual[SEQ_LENGTH+1];
	rqual[SEQ_LENGTH]='\0';

	  int i = 0;

	if (direction) {
		reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
		_tmpQual = rqual;
		_tmpSeq = _msf_seqList[readNumber].rseq;
	}
	else {
		_tmpQual = _msf_seqList[readNumber].qual;
		_tmpSeq = _msf_seqList[readNumber].seq;
	}

	int j = 0;
	int genLoc = 0;
	int *locs = (int *) l1;

	for (j = 0; j < s1; j++) {
		genLoc = locs[j];
		int af_offset[4];
		int err = -1;

		if (readFormat == 1) {
      if (key_input[index].key_number == 0) {
				if (genLoc - 1			< _msf_refGenBeg
				 || genLoc - 1 + 30 + 8 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10	- 7	< _msf_refGenBeg
				 || genLoc - 1 + 25 - 1 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 10 - 6] == readId ) {
					continue;
				}
			}
		}
		else {
      if (key_input[index].key_number == 0) {
				if (genLoc - 1 < _msf_refGenBeg
				 || genLoc - 1 + 30 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc] == readId ) {
					continue;
				}
			}
			else if (key_input[index].key_number == 1) {
				if (genLoc - 1 - 10	< _msf_refGenBeg
				 || genLoc - 1 + 20 > _msf_refGenEnd
				 || _msf_verifiedLocs[genLoc - 10 - 1] == readId ) {
					continue;
				}
			}
		}

		err = verifySingleEndSeqCG_backward(&genLoc, af_offset, _tmpSeq, key_input, index, num_base, num_key, num_ver, num_ext);

		if (err <= errThreshold && err >= 0 && _msf_verifiedLocs[genLoc] != readId) {
			generateAlignmentMatrxCG_backward(genLoc, af_offset, _tmpSeq, err, matrix);
			generateSNPSAM(matrix, strlen(matrix), editString);
			if (readFormat == 1)
				sprintf(cigar, "10M%dN10M%dN10M%dS5M", af_offset[0], af_offset[1], -af_offset[2]);
			else
				sprintf(cigar, "10M%dS10M%dS10M", -af_offset[0], -af_offset[1]);

		}
		else {
			err = -1;
		}

		if(err != -1 && !bestMode) {
			mappingCnt++;
			int offset_range = 3;
			for(i = -offset_range ; i <= offset_range ; i++) {
				if(genLoc + i >= _msf_refGenBeg && genLoc + i <= _msf_refGenEnd) {
					_msf_verifiedLocs[genLoc + i] = readId;
				}
			}

			_msf_seqList[readNumber].hits[0]++;
			_msf_output.QNAME		= _msf_seqList[readNumber].name;
			_msf_output.FLAG		= 16 * direction;
			_msf_output.RNAME		= _msf_refGenName;
			_msf_output.POS			= genLoc + _msf_refGenOffset;
			_msf_output.MAPQ		= 255;
			_msf_output.CIGAR		= cigar;
			_msf_output.MRNAME		= "*";
			_msf_output.MPOS		= 0;
			_msf_output.ISIZE		= 0;
			_msf_output.SEQ			= _tmpSeq;
			_msf_output.QUAL		= _tmpQual;
			_msf_output.optSize		= 2;
			_msf_output.optFields	= _msf_optionalFields;
			_msf_optionalFields[0].tag = "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = err;
			_msf_optionalFields[1].tag = "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = editString;

			output(_msf_output);

			if (_msf_seqList[readNumber].hits[0] == 1) {
				mappedSeqCnt++;
			}
			if ( maxHits == 0 ) {
				_msf_seqList[readNumber].hits[0] = 2;
			}

			if ( maxHits!=0 && _msf_seqList[readNumber].hits[0] == maxHits) {
				completedSeqCnt++;
				break;
			}
		}
	}
}


void resetFAST(unsigned int seqListSize) {

	freeMem(_msf_op, SEQ_LENGTH);

	if (pairedEndMode) {
		freeMem(_msf_optionalFields, 8 * sizeof(OPT_FIELDS));
	}
	else {
		freeMem(_msf_optionalFields, 2 * sizeof(OPT_FIELDS));
	}
	_msf_optionalFields = NULL;

	freeMem(_msf_samplingLocsEnds, 1);
	_msf_samplingLocsEnds = NULL;

	freeMem(_msf_sort_seqList, _msf_seqListSize * sizeof(Pair));
	freeMem(_msf_map_sort_seqList, _msf_seqListSize * sizeof(int));

	freeMem(_msf_oeaMapping, _msf_seqListSize * sizeof(int));
	freeMem(_msf_discordantMapping,  _msf_seqListSize * sizeof(int));

//	freeMem(_msf_refGenName, 4 * SEQ_LENGTH);
//printf("-- free  _msf_refGenName (%d)\n", 4 * SEQ_LENGTH);
//	_msf_refGenName = NULL;

	  if (_msf_verifiedLocs != NULL) {
			freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength + 1));
			_msf_verifiedLocs = NULL;
		}


	if (pairedEndMode) {
		freeMem(_msf_mappingInfo, seqListSize * sizeof(MappingInfo));
		freeMem(_msf_seqHits, _msf_seqListSize * sizeof(int));
		_msf_seqHits = NULL;
		freeMem(_msf_readHasConcordantMapping, _msf_seqListSize / 2 * sizeof(char));
		_msf_refGenOffset = 0;
	}
}
