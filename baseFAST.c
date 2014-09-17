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
#include "Common.h"
#include "CommandLineParser.h"
#include "Reads.h"
#include "Output.h"
#include "HashTable.h"
#include "RefGenome.h"
#include "SirFAST.h"

char *versionNumber = "2.5";
unsigned char seqFastq;
double startTime = 0;

void summaryHead();
void summaryTail();
void summaryBody();
void updateTime (double * time, double * accTime);
void updateMem (double * maxMem);


int main(int argc, char *argv[]) {
	updateTime(NULL, NULL);

	double timeLoadHash = 0.0;
	double timeLoadRead = 0.0;
	double timeMapping = 0.0;
	double timePostProc = 0.0;
	double accTimeLoadHash = 0.0;
	double accTimeLoadRead = 0.0;
	double accTimeMapping = 0.0;

	// Parsing Input Mode
	if (!parseCommandLine(argc, argv))
		return 1;
	configHashTable();

	// Indexing Mode
	if (indexingMode) {
		configHashTable();
		generateHashTable(fileName[0], fileName[1]);
		return 1;
	}

	// Searching Mode
	Read *seqList;
	unsigned int seqListSize;
	double tmpTime;;
	double maxMem = 0;

	// Loading Sequences & Sampling Locations
	int totalReads = countAllReads(seqFile1, seqFile2, seqCompressed, pairedEndMode);
	int listSize = 0;
	int accListSize = 0;
	if (totalReads <= maxInputRead)
		listSize = totalReads;
	else
		listSize = maxInputRead;

	if (pairedEndMode) {
		minPairEndedDistance = minPairEndedDistance - SEQ_LENGTH + 2;
		maxPairEndedDistance = maxPairEndedDistance - SEQ_LENGTH + 2;
	}
	if (pairedEndDiscordantMode) {
		maxPairEndedDiscordantDistance = maxPairEndedDiscordantDistance - SEQ_LENGTH + 2;
		minPairEndedDiscordantDistance = minPairEndedDiscordantDistance - SEQ_LENGTH + 2;
	}

	char outputFileName[FILE_NAME_LENGTH];
	sprintf(outputFileName, "%s%s",mappingOutputPath , mappingOutput);
	initOutput(outputFileName, outCompressed);

	if (!initLoadingHashTable(fileName[1]))
		return 1;

	fpos_t position_ih;
	fpos_t position_ig;

	updateTime(&timeLoadHash, NULL);

	if (!pairedEndMode) { // SE Mode

		FILE * _fp_unmap = initUnmapped(unmappedOutput);

		initLoadingRefGenome(fileName[0]);
		getPosHashTable(&position_ih);
		getPosRefGenome(&position_ig);

		while (readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize, listSize, accListSize)) {
			updateTime(&timeLoadRead, NULL);

			summaryHead();
			accListSize += seqListSize;

			while(loadHashTable(&tmpTime, errThreshold)) {
				initFASTCG(seqList, seqListSize, accListSize);
				updateTime(&timeLoadHash, NULL);

				mapAllSingleEndSeqCG();
				updateTime(&timeMapping, NULL);
				updateMem(&maxMem);

				summaryBody(getRefGenomeName(), timeLoadHash + timeLoadRead, timeMapping, maxMem);
				maxMem = 0;
				updateTime(&timeLoadHash, &accTimeLoadHash);
				updateTime(&timeLoadRead, &accTimeLoadRead);
				updateTime(&timeMapping, &accTimeMapping);
			}

			performUnmapped(_fp_unmap);

			summaryTail();
			setPosHashTable(&position_ih);
			setPosRefGenome(&position_ig);

			resetFAST(seqListSize);
			freeReads(seqList, seqListSize);
		}
		finalizeUnmapped(_fp_unmap);
		finalizeLoadingHashTable();
		finalizeFAST();
	}
	else { // PE Mode
		int unmappedCnt = 0;
		int unmappedCnt_slice = 0;

		FILE * _fp_unmap = initUnmapped(unmappedOutput);	// output for unmapped reads
		FILE * _fp_oea = initOEAReads(outputFileName);		// output for oea reads
		FILE * _fp_divet = initPairedEndDiscPP();					// output for divet

		initLoadingRefGenome(fileName[0]);
		getPosHashTable(&position_ih);
		getPosRefGenome(&position_ig);

		while (readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize, listSize, accListSize)) {
			updateTime(&timeLoadRead, NULL);
			initBestMapping(seqListSize);
			summaryHead();
			accListSize += seqListSize;

			while(loadHashTable(&tmpTime, errThreshold)) {
				initFASTCG(seqList, seqListSize, accListSize);
				updateTime(&timeLoadHash, NULL);
				mapPairedEndSeqCG();
				updateTime(&timeMapping, NULL);
				unmappedCnt = outputPairedEnd(unmappedCnt_slice);	// output for unmapped reads
				updateTime(&timePostProc, NULL);
				updateMem(&maxMem);
        summaryBody(getRefGenomeName(), timeLoadHash + timeLoadRead, timeMapping, maxMem);
				maxMem = 0;
				updateTime(&timeLoadHash, &accTimeLoadHash);
				updateTime(&timeLoadRead, &accTimeLoadRead);
				updateTime(&timeMapping, &accTimeMapping);
			}
			unmappedCnt_slice = unmappedCnt;

			performUnmapped(_fp_unmap);						// output for unmapped reads
			performOEAReads(_fp_oea);							// output for oea
			if (pairedEndDiscordantMode)
				performPairedEndDiscPP(_fp_divet);	// output for divet
			summaryTail();
			setPosHashTable(&position_ih);
			setPosRefGenome(&position_ig);
			performBestMapping(seqListSize);
			finalizeBestMapping(seqListSize);
			resetFAST(seqListSize);
			freeReads(seqList, seqListSize);
			updateTime(&timePostProc, NULL);
		}

		finalizeUnmapped(_fp_unmap);
		finalizeOEAReads(_fp_oea);
		if (pairedEndDiscordantMode)
			finalizePairedEndDiscPP(_fp_divet);
		outputAllTransChromosomal(transChromosomal);	// not implemented
		finalizeLoadingHashTable();
		finalizeFAST();
	}

	finalizeOutput();
	updateTime(&timePostProc, NULL);

  fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:", accTimeLoadHash + accTimeLoadRead, accTimeMapping);
	fprintf(stdout, "Post Processing Time: %18.2f \n", timePostProc);
	fprintf(stdout, "%-30s%10.2f\n","Total Time:", accTimeLoadHash + accTimeLoadRead + accTimeMapping + timePostProc);
	fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", accListSize);
	fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
	fprintf(stdout, "%-30s%10.0f\n\n","Avg No. of locations verified:", ceil((float)verificationCnt/seqListSize));

	int cof = (pairedEndMode) ? 2 : 1;
	if (progressRep && maxHits != 0) {
    int frequency[maxHits + 1];
    int i;
    for (i = 0 ; i <= maxHits; i++)
        frequency[i] = 0;
    for (i = 0; i < seqListSize; i++)
			frequency[(int)(*(seqList[i*cof].hits))]++;
		frequency[maxHits] = completedSeqCnt;
		for (i = 0 ; i <= maxHits; i++)
			fprintf(stdout, "%-30s%10d%10d%10.2f%%\n","Reads Mapped to ", i, frequency[i], 100*(float)frequency[i]/(float)seqListSize);
	}
  return 0;
}

void summaryHead(){
  fprintf(stdout, "-------------------------------------------------------");
	fprintf(stdout, "------------------------------------------------------\n");
	fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s | %15s |\n", "Seq. Name",
    "Loading Time", "Mapping Time", "Memory Usage(M)", "Total Mappings", "Mapped reads");
	fprintf(stdout, "-------------------------------------------------------");
	fprintf(stdout, "------------------------------------------------------\n");
	fflush(stdout);
}

void summaryTail() {
	fprintf(stdout, "-------------------------------------------------------");
	fprintf(stdout, "------------------------------------------------------\n");
	fflush(stdout);
}

void summaryBody (char* curGen, double loadingTime, double mappingTime, double maxMem) {
	fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld | %15lld |\n",
						curGen, loadingTime, mappingTime, maxMem, mappingCnt, mappedSeqCnt);
	fflush(stdout);
}

void updateTime (double * time, double * accTime) {
  if (time != NULL)
    *time += getTime() - startTime;

  if (accTime != NULL) {
    *accTime += *time;
    *time = 0;
  }

  startTime = getTime();
}

void updateMem (double * maxMem) {
  if (*maxMem < getMemUsage()) {
    *maxMem = getMemUsage();
  }
}
