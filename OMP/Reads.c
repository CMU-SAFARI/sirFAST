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
#include <ctype.h>
#include <zlib.h>
#include "Common.h"
#include "Reads.h"

#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))

FILE *_r_fp1;
FILE *_r_fp2;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
int _r_seqCnt;
int *_r_samplingLocs;

/**********************************************/



/**********************************************/

char *(*readFirstSeq)(char *);

char *(*readSecondSeq)(char *);

/**********************************************/
char *readFirstSeqTXT( char *seq )
{
  return fgets(seq, SEQ_MAX_LENGTH, _r_fp1);
}

/**********************************************/
char *readSecondSeqTXT( char *seq )
{
  return fgets(seq, SEQ_MAX_LENGTH, _r_fp2);
}
/**********************************************/
char *readFirstSeqGZ( char *seq )
{
  return gzgets(_r_gzfp1, seq, SEQ_MAX_LENGTH);
}

/**********************************************/
char *readSecondSeqGZ( char *seq )
{
  return gzgets(_r_gzfp2, seq, SEQ_MAX_LENGTH);
}

/**********************************************/
int toCompareRead(const void * elem1, const void * elem2)
{
  return strcmp(((Read *)elem1)->seq, ((Read *)elem2)->seq);	
}

/**********************************************/
int readAllReads(char *fileName1, char *fileName2, int compressed,
		 unsigned char *fastq, unsigned char pairedEnd, Read **seqList,
		 unsigned int *seqListSize, unsigned int ListSize, unsigned int AccListSize) {

	double startTime=getTime();
	
	char * seq;
	char * qual;
	char seq1[SEQ_MAX_LENGTH];
	char rseq1[SEQ_MAX_LENGTH];
	char qual1[SEQ_MAX_LENGTH];
	char seq2[SEQ_MAX_LENGTH];
	char rseq2[SEQ_MAX_LENGTH];
	char qual2[SEQ_MAX_LENGTH];
	char dummy[SEQ_MAX_LENGTH];
	int discarded = 0;
	int seqCnt = 0;
	Read *list = NULL;
	int nCnt1;
	int nCnt2;

	list = getMem(sizeof(Read)*ListSize, "list @readAllReads()");

	while(ListSize > seqCnt && readFirstSeq(dummy)) {
		int i = 0;
		int _mtmp = 36;

		if(dummy[0] == '#' || dummy[0] == '>' || dummy[0] == ' ' || 
			dummy[0] == '\r' || dummy[0] == '\n')
			continue;

		strtok(dummy, "\t ");
		seq = strtok(NULL, "\t ");
		qual = strtok(NULL, "\t ");

		for(i = 0; i < _mtmp - 1; i++) {
			seq1[i] = toupper(seq[i]);
			qual1[i] = qual[i];
		}

		for(i = 0; i < _mtmp - 1; i++) {
			seq2[i] = toupper(seq[i + _mtmp - 1]);
			qual2[i] = qual[i + _mtmp - 1];
		}

		seq1[_mtmp - 1] = seq2[_mtmp - 1] = qual1[_mtmp - 1] = qual2[_mtmp - 1] = '\0';
		
		nCnt1 = 0; nCnt2 = 0;
		for (i=0; i<_mtmp; i++)
		  {
		    if (seq1[i] == 'N')		      
			nCnt1++;
		    if (seq2[i] == 'N')		      
			nCnt2++;		    
		  }
		
		if (nCnt1 > errThreshold || nCnt2 > errThreshold)
		  {
		    discarded += 2;
		    continue;
		  }
	

		if (errThreshold == 255) {
			if (cropSize > 0) {
				errThreshold = (int) ceil(cropSize * 0.04);
				fprintf(stdout, "Sequence length: %d bp. Error threshold is set to %d bp.\n", 
					cropSize, errThreshold);
				}
			else {
				errThreshold = (int) ceil((strlen(seq1)) * 0.04);
				fprintf(stdout, "Sequence length: %d bp. Error threshold is set to %d bp.\n", 
					((int)strlen(seq1)), errThreshold);
			}
				fprintf(stdout, "You can override this value using the -e parameter.\n");
		}

		list[seqCnt].hits = getMem (1 + 3 * _mtmp + 3 + _mtmp, "list.hits @readAllReads()");
		list[seqCnt].seq = list[seqCnt].hits + 1;
		list[seqCnt].rseq = list[seqCnt].seq + _mtmp + 1;
		list[seqCnt].qual = list[seqCnt].rseq + _mtmp + 1;
		list[seqCnt].name = list[seqCnt].qual + _mtmp + 1;
		list[seqCnt].hashValue = getMem(sizeof(short) * _mtmp, "list.hashValue @readAllReads()");
		list[seqCnt].rhashValue = getMem(sizeof(short) * _mtmp, "list.rhashValue @readAllReads()");
		list[seqCnt].readNumber = seqCnt;			
		list[seqCnt].hits[0] = 0;

		reverseComplement(seq1, rseq1, _mtmp - 1);	// DHL Modify
		rseq1[_mtmp - 1] = '\0';

		for (i=0; i<_mtmp-1; i++) {
			list[seqCnt].seq[i] = seq1[i];
			list[seqCnt].rseq[i] = rseq1[i];
			list[seqCnt].qual[i] = qual1[i];
		}

		if (!pairedEndMode)
		  sprintf(list[seqCnt].name, "%s_%d/1", mappingOutput, (seqCnt + AccListSize) / 2);
		else
		  sprintf(list[seqCnt].name, "%s_%d", mappingOutput, (seqCnt + AccListSize) / 2);

		list[seqCnt].seq[_mtmp - 1] = list[seqCnt].rseq[_mtmp - 1] = list[seqCnt].qual[_mtmp - 1]='\0';

		seqCnt++;

		list[seqCnt].hits = getMem (1 + 3 * _mtmp + 3 + _mtmp, "list.hits @readAllReads()");
		list[seqCnt].seq = list[seqCnt].hits + 1;
		list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
		list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
		list[seqCnt].name = list[seqCnt].qual + _mtmp+1;
		list[seqCnt].hashValue = getMem(sizeof(short) * _mtmp, "list.hashValue @readAllReads()");
		list[seqCnt].rhashValue = getMem(sizeof(short) * _mtmp, "list.rhashValue @readAllReads()");
		list[seqCnt].readNumber = seqCnt;				 
		list[seqCnt].hits[0] = 0;

		reverseComplement(seq2, rseq2, _mtmp - 1);	// DHL Modify
		rseq2[_mtmp - 1] = '\0';

		for (i=0; i<_mtmp; i++) {
			list[seqCnt].seq[i] = seq2[i];
			list[seqCnt].rseq[i] = rseq2[i];
			list[seqCnt].qual[i] = qual2[i];
		}

		if (!pairedEndMode)
		  sprintf(list[seqCnt].name, "%s_%d/2", mappingOutput, (seqCnt + AccListSize) / 2);
		else
		  sprintf(list[seqCnt].name, "%s_%d", mappingOutput, (seqCnt + AccListSize) / 2);

		list[seqCnt].seq[_mtmp - 1] = list[seqCnt].rseq[_mtmp - 1] = list[seqCnt].qual[_mtmp - 1]='\0';

		seqCnt++;
	}

	if (seqCnt <= 0) {
		//fprintf(stdout, "ERROR: No reads can be found for mapping\n");
		fprintf(stdout, "==== End of Input Reads    					====\n");	// DHL: read slice
		return 0;
	}

	//qsort(list, seqCnt, sizeof(Read), toCompareRead);
	adjustQual(list, seqCnt);
	*seqList = list;
	*seqListSize = seqCnt;
	_r_seq = list;
	_r_seqCnt = seqCnt;

	if (pairedEnd) 
		discarded *= 2;

	if (seqCnt > 1) {
		fprintf(stdout, "==== %d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]	====\n", 
			seqCnt, (getTime()-startTime), discarded, getMemUsage());
	}
	else {
		fprintf(stdout, "==== %d sequence is read in %0.2f. (%d discarded) [Mem:%0.2f M]	====\n", 
			seqCnt, (getTime()-startTime), discarded, getMemUsage());
	}

	return seqCnt;
}

/********************************************************************************************/
void freeReads(Read *seqList, unsigned int seqListSize) {
	int seqCnt = 0;
	int _mtmp = 36;
//	char seq1[SEQ_MAX_LENGTH];

	for (seqCnt = 0; seqCnt < seqListSize; seqCnt++) {
		//freeMem(seqList[seqCnt].hits, 1 + 3 * _mtmp + 3 + strlen(seq1) + 1, "_seqList.hits @freeReads()");
		freeMem(seqList[seqCnt].hits, 1 + 3 * _mtmp + 3 + _mtmp, "_seqList.hits @freeReads()");
		freeMem(seqList[seqCnt].hashValue, sizeof(short) * _mtmp, "_seqList.hashValue @freeReads()");
		freeMem(seqList[seqCnt].rhashValue, sizeof(short) * _mtmp, "_seqList.rhashValue @freeReads()");
	}
	freeMem(seqList, sizeof(Read) * seqListSize, "seqList @freeReads()");
	return;
}


/********************************************************************************************/
int countAllReads(char *fileName1, char *fileName2, int compressed, 
		unsigned char pairedEnd) {


	char dummy[SEQ_MAX_LENGTH];
	int maxCnt = 0;
	
	if (!compressed) {
		_r_fp1 = fileOpen( fileName1, "r");
		if (_r_fp1 == NULL)
			return 0;

		if ( pairedEnd && fileName2 != NULL ) {
			_r_fp2 = fileOpen ( fileName2, "r" );
			if (_r_fp2 == NULL)
				return 0;
		}
		else {
			_r_fp2 = _r_fp1;
		}

		readFirstSeq = &readFirstSeqTXT;
		readSecondSeq = &readSecondSeqTXT;
	}
	else {
		_r_gzfp1 = fileOpenGZ (fileName1, "r");
		if (_r_gzfp1 == NULL)
			return 0;

		if ( pairedEnd && fileName2 != NULL ) {
			_r_gzfp2 = fileOpenGZ ( fileName2, "r" );
			if (_r_gzfp2 == NULL)
				return 0;
		}
		else {
			_r_gzfp2 = _r_gzfp1;
		}

		readFirstSeq = &readFirstSeqGZ;
		readSecondSeq = &readSecondSeqGZ;
	}

	// Counting the number of lines in the file
	while (readFirstSeq(dummy)) { 
		if(dummy[0] != '#' && dummy[0]!='>' && dummy[0] != ' ' && 
			dummy[0] != '\r' && dummy[0] != '\n')
			maxCnt++;
	}

	if (!compressed)
		rewind(_r_fp1);
	else
		gzrewind(_r_gzfp1);

	// Return the Maximum # of sequences
	return maxCnt * 2;
}

/********************************************************************************************/
void closingReads (char *fileName2, unsigned char pairedEnd, int compressed) {
	if (!compressed) {
		fclose(_r_fp1);
		if (pairedEnd && fileName2 != NULL) 
			fclose(_r_fp2);
	}
	return;

	gzclose(_r_gzfp1);
	if ( pairedEnd && fileName2 != NULL)
		gzclose(_r_gzfp2);
	return;
}


/**********************************************/
void loadSamplingLocations(int **samplingLocs, int * samplingLocsSize)
{
  int i;
  int samLocsSize = errThreshold + 1;
  int *samLocs = getMem(sizeof(int)*samLocsSize, "samLocs @loadSamplingLocations()");

  for (i=0; i<samLocsSize; i++)
    {
      samLocs[i] = (SEQ_LENGTH / samLocsSize) *i;
      if ( samLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
	samLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
    }

  // Outputing the sampling locations

  /*

    int j;
    for (i=0; i<SEQ_LENGTH; i++)
    {
    fprintf(stdout, "-");
    }
    fprintf(stdout, "\n");

    for ( i=0; i<samLocsSize; i++ )
    {
    for ( j=0; j<samLocs[i]; j++ )
    {
    fprintf(stdout," ");
    }
    for (j=0; j<WINDOW_SIZE; j++)
    {
    fprintf(stdout,"+");
    }
    fprintf(stdout, "\n");
    fflush(stdout);
    }
	

    for ( i=0; i<SEQ_LENGTH; i++ )
    {
    fprintf(stdout, "-");
    }
    fprintf(stdout, "\n");

  */

  *samplingLocs = samLocs;
  *samplingLocsSize = samLocsSize;
  _r_samplingLocs = samLocs;
}

void adjustQual(Read *list, int seqCnt){
  /* This function will automatically determine the phred_offset and readjust quality values if needed */
  int i,j,q, offset=64;
  int len = strlen(list[0].qual);
  
  for (i=0; i<10000 && i<seqCnt; i++){
    for (j=0;j<len;j++){
      q = (int) list[i].qual[j] - offset;
      if (q < 0){
	offset = 33;
	break;
      }
    }
    if (offset == 33)
      break;
  }
  
  if (offset == 64){
    fprintf(stdout, "[Quality Warning] Phred offset is 64. Readjusting to 33.\n");
    fflush(stdout);
    for (i=0;i<seqCnt;i++){
      for (j=0;j<len;j++){
	list[i].qual[j] -= 31;
      }
    }
  }
}

FILE * initUnmapped(char *fileName) {
	FILE * _fp_unmap = NULL;
	if (fileName != NULL)
		_fp_unmap = fileOpen(fileName, "w");
	return _fp_unmap;
}

void finalizeUnmapped(FILE * _fp_unmap) {
	fclose(_fp_unmap);
}

void operateUnmapped(FILE * fp1) {

	if (pairedEndMode)
		_r_seqCnt /=2;

	int i = 0;

	for (i = 0; i < _r_seqCnt; i++) {
		if (pairedEndMode && _r_seq[2*i].hits[0] == 0 && _r_seq[2*i+1].hits[0] == 0	&&	strcmp(_r_seq[2*i].qual, "*") != 0)
			fprintf(fp1,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
		else if (pairedEndMode && _r_seq[2*i].hits[0] == 0 && _r_seq[2*i+1].hits[0] == 0)
			fprintf(fp1, ">%s/1\n%s\n>%s/2\n%shits=%d\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[2*i+1].hits[0]);
		else if (!pairedEndMode && _r_seq[i].hits[0] == 0 && strcmp(_r_seq[i].qual, "*")!=0)
			fprintf(fp1,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
		else if (!pairedEndMode && _r_seq[i].hits[0] == 0)
			fprintf(fp1,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
	}
	if (pairedEndMode)
		_r_seqCnt *= 2;

//	for (i = 0; i < _r_seqCnt; i++)
//		freeMem(_r_seq[i].hits,0);
//	freeMem(_r_seq,0);
//	freeMem(_r_samplingLocs,0);	// Modified by DHL
}
