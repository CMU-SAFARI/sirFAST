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
 *	 list of conditions and the following disclaimer in the documentation and/or other
 *	 materials provided with the distribution.
 * - Neither the names of the University of Washington, Simon Fraser University, 
 *	 nor the names of its contributors may be
 *	 used to endorse or promote products derived from this software without specific
 *	 prior written permission.
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

char 				*versionNumber = "2.5";			// Current Version
unsigned char		seqFastq;

int main(int argc, char *argv[]) {

	if (!parseCommandLine(argc, argv))
	return 1;

	configHashTable();

	/****************************************************
	 * INDEXING
	 ***************************************************/
	if (indexingMode) {
		configHashTable();
		generateHashTable(fileName[0], fileName[1]);
		return 1;
	}
	/****************************************************
	 * SEARCHING
	 ***************************************************/
	Read *seqList;
	unsigned int seqListSize;
	int fc;
	int samplingLocsSize;
	int *samplingLocs;
	double totalLoadingTime = 0;
	double totalMappingTime = 0;
	double startTime;
	double loadingTime;
	double mappingTime;
	double lstartTime;
	double ppTime = 0.0;
	double tmpTime;;
	char *prevGen = getMem(CONTIG_NAME_SIZE, "prevGen @baseFAST()");
	prevGen[0] = '\0';
	char *curGen;
	int	flag;
	double maxMem=0;
	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];
	char outputFileName[FILE_NAME_LENGTH];

	int TotalReads = 0;
	int readsCnt = 0;
	int ListSize = maxInputRead;
	int AccReads = 0;
	int pre_unmappedCnt = 0;
	int pre_unmappedCnt_slice = 0;
	fpos_t position_ih;
	fpos_t position_rg;

	// Loading Sequences & Sampling Locations
	startTime = getTime();																		// LOAD TIME (READ)
	TotalReads = countAllReads(seqFile1, seqFile2, seqCompressed, pairedEndMode);				// LOAD TIME (READ)
	if (TotalReads <= maxInputRead)
		ListSize = TotalReads;
	fprintf(stdout, "\n==== Total number of sequences: %d				====\n", TotalReads);	// LOAD TIME (READ)
	readsCnt = readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode,		// LOAD TIME (READ)
		&seqList, &seqListSize, ListSize, AccReads);											// LOAD TIME (READ)
	AccReads += seqListSize;																	// LOAD TIME (READ)
	if (!readsCnt)																				// LOAD TIME (READ)
		return 1;																				// LOAD TIME (READ)
	loadSamplingLocations(&samplingLocs, &samplingLocsSize);									// LOAD TIME (READ)
	totalLoadingTime += getTime() - startTime;													// LOAD TIME (READ)

	if (pairedEndMode) {
		minPairEndedDistance = minPairEndedDistance - SEQ_LENGTH + 2;
		maxPairEndedDistance = maxPairEndedDistance - SEQ_LENGTH + 2;
		if (pairedEndDiscordantMode) {
			maxPairEndedDiscordantDistance = maxPairEndedDiscordantDistance - SEQ_LENGTH + 2;
			minPairEndedDiscordantDistance = minPairEndedDiscordantDistance - SEQ_LENGTH + 2;
		}
		sprintf(fname1, "__%s__1", mappingOutput);
		sprintf(fname2, "__%s__2", mappingOutput);
		sprintf(fname3, "__%s__disc", mappingOutput);
		sprintf(fname4, "__%s__oea1", mappingOutput);
		sprintf(fname5, "__%s__oea2", mappingOutput);
		unlink(fname1);
		unlink(fname2);
		unlink(fname3);
		unlink(fname4);
		unlink(fname5);
	}

	sprintf(outputFileName, "%s%s", mappingOutputPath, mappingOutput);
	initOutput(outputFileName, outCompressed);
	fprintf(stdout, "Number of threads: %d\n", number_of_threads);

	fprintf(stdout, "-------------------------------------------------------");
	fprintf(stdout, "------------------------------------------------------\n");
	fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s | %15s |\n", "Seq. Name", 
		"Loading Time", "Mapping Time", "Memory Usage(M)", "Total Mappings", "Mapped reads");
	fprintf(stdout, "-------------------------------------------------------");
	fprintf(stdout, "------------------------------------------------------\n");
	fflush(stdout);

	// Single End Mode
	if (!pairedEndMode) {
		if(bestMode) {
			initBestMapping(readsCnt);
			//resetFAST(readsCnt);
		}
		if (!initLoadingHashTable(fileName[1]))
			return 1;
		
		sprintf(outputFileName, "%s%s", mappingOutputPath, mappingOutput);		
		initFASTCG(seqList, seqListSize, fileName[0], AccReads, 1);

		FILE * _fp_unmap = NULL;
		_fp_unmap = initUnmapped(unmappedOutput);

		getPosHashTable(&position_ih);

		mappingTime = 0;
		loadingTime = 0;
		prevGen[0] = '\0';
		flag = 1;

		do {
			flag = loadHashTable ( &tmpTime, errThreshold);				// Reading a fragment
			curGen = getRefGenomeName();

			if (flag && prevGen[0]== '\0') 
				sprintf(prevGen, "%s", curGen);

			if ( !flag || strcmp(prevGen, curGen)!=0){

				fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld | %15lld |\n", 
					prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
				fflush(stdout);

				totalMappingTime += mappingTime;
				totalLoadingTime += loadingTime;
				loadingTime = 0;
				mappingTime = 0;
				maxMem = 0;

				if (!flag) {

					operateUnmapped(_fp_unmap);
					if (bestMode) {
						finalizeBestSingleMapping();
					}

					resetFAST(seqListSize);
					freeReads(seqList, seqListSize);

					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n\n");

					startTime = getTime();														// LOAD TIME (READ)
					if (TotalReads - AccReads < ListSize) ListSize = TotalReads - AccReads;		// LOAD TIME (READ)
					readsCnt = readAllReads(seqFile1, seqFile2, seqCompressed,					// LOAD TIME (READ)
						&seqFastq, pairedEndMode, &seqList, &seqListSize, ListSize, AccReads);	// LOAD TIME (READ)
					if (!readsCnt) break;														// LOAD TIME (READ)
					totalLoadingTime += getTime() - startTime;									// LOAD TIME (READ)
	
					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n");
					fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s | %15s |\n", "Seq. Name", 
						"Loading Time", "Mapping Time", "Memory Usage(M)", "Total Mappings", "Mapped reads");
					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n");
					AccReads += seqListSize;
				
					if (bestMode) {	
						initBestMapping(readsCnt);
					}

					//updateSeqList(seqListSize, seqList);
	
					setPosHashTable(&position_ih);
	
					flag = loadHashTable(&tmpTime, errThreshold);	// Load TIME
					loadingTime += tmpTime;							// Load TIME
	
					curGen = getRefGenomeName();
					if (flag && prevGen[0] == '\0') 
						sprintf(prevGen, "%s", curGen);
				}
			}
			sprintf(prevGen, "%s", curGen);
			loadingTime += tmpTime;
			initFASTCG(seqList, seqListSize, fileName[0], AccReads, 0);
			lstartTime = getTime();
			mapAllSingleEndSeqCG();	
			mappingTime += getTime() - lstartTime;
			if (maxMem < getMemUsage()) maxMem = getMemUsage();					 
		} while (flag);

		finalizeUnmapped(_fp_unmap);

//		if(bestMode)
//			finalizeBestSingleMapping();
//		finalizeFAST();
		finalizeLoadingHashTable();
	}
	// Paired End Mode
	else {
		if (pairedEndMode) {
			initBestMapping(readsCnt);
	//		resetFAST(readsCnt);
		}
		if (!initLoadingHashTable(fileName[1]))
			return 1;
	
		sprintf(outputFileName, "%s%s", mappingOutputPath, mappingOutput);		
		initFASTCG(seqList, seqListSize, fileName[0], AccReads, 1);
	
		FILE * _fp_oea = initOEAReads(outputFileName);
		FILE * _fp_unmap = initUnmapped(unmappedOutput);
		FILE * _fp_divet = initPairedEndDiscPP();
	
		getPosHashTable(&position_ih);
		getPosRefGenome(&position_rg);
	
		mappingTime = 0;
		loadingTime = 0;
		ppTime = 0;
		prevGen[0] = '\0';
		flag = 1;
	
		do {
			flag = loadHashTable(&tmpTime, errThreshold);	// Load TIME
			loadingTime += tmpTime;							// Load TIME
			curGen = getRefGenomeName();
			if (flag && prevGen[0] == '\0') sprintf(prevGen, "%s", curGen);
			if (!flag || strcmp(prevGen, curGen) != 0) {
				lstartTime = getTime();										// MAPPING TIME
				pre_unmappedCnt = outputPairedEnd(pre_unmappedCnt_slice);	// MAPPING TIME
				mappingTime += getTime() - lstartTime;						// MAPPING TIME
	
				fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld | %15lld |\n", 
					prevGen, loadingTime, mappingTime, maxMem, mappingCnt, mappedSeqCnt);
				fflush(stdout);

				totalMappingTime += mappingTime;
				totalLoadingTime += loadingTime;
				loadingTime = 0;
				mappingTime = 0;
				maxMem = 0;
	
				if (!flag) {
					pre_unmappedCnt_slice = pre_unmappedCnt;	
					operateOEAReads(_fp_oea);
					outputAllTransChromosomal(transChromosomal);
					finalizeBestConcordantDiscordant();
	
					if (pairedEndDiscordantMode) {
						lstartTime = getTime();				// PPTIME
						operatePairedEndDiscPP(_fp_divet);	// PPTIME
						ppTime += getTime() - lstartTime;	// PPTIME
					}
	
					operateUnmapped(_fp_unmap);
	
					resetFAST(seqListSize);
					freeReads(seqList, seqListSize);
					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n\n");
	
					startTime = getTime();														// LOAD TIME (READ)
					if (TotalReads - AccReads < ListSize) ListSize = TotalReads - AccReads;		// LOAD TIME (READ)
					readsCnt = readAllReads(seqFile1, seqFile2, seqCompressed,					// LOAD TIME (READ)
						&seqFastq, pairedEndMode, &seqList, &seqListSize, ListSize, AccReads);	// LOAD TIME (READ)
					if (!readsCnt) break;														// LOAD TIME (READ)
					totalLoadingTime += getTime() - startTime;									// LOAD TIME (READ)
	
					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n");
					fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s | %15s |\n", "Seq. Name", 
						"Loading Time", "Mapping Time", "Memory Usage(M)", "Total Mappings", "Mapped reads");
					fprintf(stdout, "-------------------------------------------------------");
					fprintf(stdout, "------------------------------------------------------\n");
					AccReads += seqListSize;
	
					initBestMapping(seqListSize);
					//updateSeqList(seqListSize, seqList);
	
					setPosHashTable(&position_ih);
					setPosRefGenome(&position_rg);	
	
					flag = loadHashTable(&tmpTime, errThreshold);	// Load TIME
					loadingTime += tmpTime;							// Load TIME
	
					curGen = getRefGenomeName();
					if (flag && prevGen[0] == '\0') 
						sprintf(prevGen, "%s", curGen);
				}
			}
	
			sprintf(prevGen, "%s", curGen);
			initFASTCG(seqList, seqListSize, fileName[0], AccReads, 0);
			mapPairedEndSeqCG(seqList, seqListSize, AccReads);

			if (maxMem < getMemUsage()) 
				maxMem = getMemUsage();
		} while (flag);


		if (pairedEndDiscordantMode) {
			finalizePairedEndDiscPP(_fp_divet);
		}

		finalizeOEAReads(_fp_oea);	
		finalizeUnmapped(_fp_unmap);
		finalizeLoadingHashTable();
		finalizeFAST();
		closingReads(seqFile2, pairedEndMode, seqCompressed);
	}

	finalizeOutput(outputFileName);

	fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:", totalLoadingTime, totalMappingTime);
	if (pairedEndDiscordantMode)
		fprintf(stdout, "Post Processing Time: %18.2f \n", ppTime);
	fprintf(stdout, "%-30s%10.2f\n", "Total Time:", totalMappingTime+totalLoadingTime);
	fprintf(stdout, "%-30s%10d\n", "Total No. of Reads:", AccReads);
	fprintf(stdout, "%-30s%10lld\n", "Total No. of Mappings:", mappingCnt);
	fprintf(stdout, "%-30s%10.0f\n\n", "Avg No. of locations verified:", 
		ceil((float)verificationCnt/seqListSize));

	int cof = (pairedEndMode)?2:1;

	if (progressRep && maxHits != 0) {
		int frequency[maxHits+1];
		int i;
		for (i = 0; i <= maxHits; i++) {
			frequency[i] = 0;
		}

		for (fc = 0; fc < seqListSize; fc++) {
			frequency[(int)(*(seqList[fc*cof].hits))]++;
		}

		frequency[maxHits] = completedSeqCnt;

		for (i = 0; i <= maxHits; i++) {
			fprintf(stdout, "%-30s%10d%10d%10.2f%%\n", "Reads Mapped to ", 
				i, frequency[i], 100 * (float)frequency[i] / (float)seqListSize);
		}
	}

	freeMem(prevGen, CONTIG_NAME_SIZE, "prevGen @baseFAST()");
	return 0;
}
