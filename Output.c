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
#include <zlib.h>
#include <string.h>
#include "Common.h"
#include "Output.h"

FILE			*_out_fp;
gzFile			_out_gzfp;

char buffer[300000];
int bufferSize = 0;

void finalizeGZOutput() {
  gzclose(_out_gzfp);
}

void finalizeTXOutput() {
  fclose(_out_fp);
}


void gzOutputQ(SAM map) {
	int i = 0;
	int index = 0;
	int gap1 = 0;
	int gap2 = 0;
	int gap3 = 0;

	char gc_cg[SEQ_LENGTH];
	char md_cg[SEQ_LENGTH];
	char seq_cg[SEQ_LENGTH];
	char gs_cg[SEQ_LENGTH];
	char qual_cg[SEQ_LENGTH];
	char gq_cg[SEQ_LENGTH];

	if(sscanf(map.CIGAR, "5M%dS10M%dN10M%dN10M", &gap1, &gap2, &gap3) == 3) {
		if(gap2 == 0)
			sprintf(md_cg, "%dM%dN10M", 25-gap1, gap3);
		else
			sprintf(md_cg, "%dM%dN10M%dN10M", 15-gap1,gap2, gap3);

		sprintf(gc_cg, "%dS%dG%dS", 5-gap1, gap1, 30-gap1);
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i < 5 || i >= 5+gap1) {
				seq_cg[index] = map.SEQ[i];
				qual_cg[index] = map.QUAL[i];
				index++;
			}
		}
		seq_cg[index]='\0';
		qual_cg[index]='\0';
		index = 0;
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i >= 5-gap1 && i < 5+gap1) {
				gs_cg[index] = map.SEQ[i];
				gq_cg[index] = map.QUAL[i];
				index++;
			}
		}
		gs_cg[index]='\0';
		gq_cg[index]='\0';
	}
	else {
		sscanf(map.CIGAR, "10M%dN10M%dN10M%dS5M", &gap1, &gap2, &gap3);
		if(gap2 == 0)
			sprintf(md_cg, "10M%dN%dM", gap1, 25-gap3);
		else
			sprintf(md_cg, "10M%dN10M%dN%dM", gap1,gap2, 15-gap3);
		sprintf(gc_cg, "%dS%dG%dS", 30-gap3, gap3, 5-gap3);
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i < 30 || i >= 30+gap3) {
				seq_cg[index] = map.SEQ[i];
				qual_cg[index] = map.QUAL[i];
				index++;
			}
		}
		seq_cg[index]='\0';
		qual_cg[index]='\0';
		index = 0;
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i >= 30-gap3 && i < 30+gap3) {
				gs_cg[index] = map.SEQ[i];
				gq_cg[index] = map.QUAL[i];
				index++;
			}
		}
		gs_cg[index]='\0';
		gq_cg[index]='\0';
	}

	gzprintf(_out_gzfp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			map.QNAME, map.FLAG, map.RNAME, map.POS+1, map.MAPQ, md_cg, map.MRNAME, map.MPOS, map.ISIZE, seq_cg, qual_cg);

	for ( i = 0; i < map.optSize; i++) {
		switch (map.optFields[i].type) {
			case 'A':
				gzprintf(_out_gzfp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
				break;
			case 'i':
				gzprintf(_out_gzfp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
				break;
			case 'f':
				gzprintf(_out_gzfp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
				break;
			case 'Z':
			case 'H':
				gzprintf(_out_gzfp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
				break;
		}
	}

	gzprintf(_out_gzfp, "\tGC:Z:%s", gc_cg);
	gzprintf(_out_gzfp, "\tGS:Z:%s", gs_cg);
	gzprintf(_out_gzfp, "\tGQ:Z:%s", gq_cg);
	if (readGroup[0] != 0)
		gzprintf(_out_gzfp, "\t%s:%c:%s", "RG", 'Z', readGroup);

	gzprintf(_out_gzfp, "\n");
}

void outputSAM(FILE *fp, SAM map) {
	int i = 0;
	int index = 0;
	int gap1 = 0;
	int gap2 = 0;
	int gap3 = 0;

	char gc_cg[SEQ_LENGTH];
	char md_cg[SEQ_LENGTH];
	char seq_cg[SEQ_LENGTH];
	char gs_cg[SEQ_LENGTH];
	char qual_cg[SEQ_LENGTH];
	char gq_cg[SEQ_LENGTH];

	if(sscanf(map.CIGAR, "5M%dS10M%dN10M%dN10M", &gap1, &gap2, &gap3) == 3) {
		if(gap2 == 0)
			sprintf(md_cg, "%dM%dN10M", 25-gap1, gap3);
		else
			sprintf(md_cg, "%dM%dN10M%dN10M", 15-gap1,gap2, gap3);
		sprintf(gc_cg, "%dS%dG%dS", 5-gap1, gap1, 30-gap1);
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i < 5 || i >= 5+gap1) {
				seq_cg[index] = map.SEQ[i];
				qual_cg[index] = map.QUAL[i];
				index++;
			}
		}
		seq_cg[index]='\0';
		qual_cg[index]='\0';
		index = 0;
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i >= 5-gap1 && i < 5+gap1) {
				gs_cg[index] = map.SEQ[i];
				gq_cg[index] = map.QUAL[i];
				index++;
			}
		}
		gs_cg[index]='\0';
		gq_cg[index]='\0';
	}
	else {
		sscanf(map.CIGAR, "10M%dN10M%dN10M%dS5M", &gap1, &gap2, &gap3);
		if(gap2 == 0)
			sprintf(md_cg, "10M%dN%dM", gap1, 25-gap3);
		else
			sprintf(md_cg, "10M%dN10M%dN%dM", gap1,gap2, 15-gap3);
		sprintf(gc_cg, "%dS%dG%dS", 30-gap3, gap3, 5-gap3);
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i < 30 || i >= 30+gap3) {
				seq_cg[index] = map.SEQ[i];
				qual_cg[index] = map.QUAL[i];
				index++;
			}
		}
		seq_cg[index]='\0';
		qual_cg[index]='\0';
		index = 0;
		for(i = 0; i < SEQ_LENGTH; i++) {
			if(i >= 30-gap3 && i < 30+gap3) {
				gs_cg[index] = map.SEQ[i];
				gq_cg[index] = map.QUAL[i];
				index++;
			}
		}
		gs_cg[index]='\0';
		gq_cg[index]='\0';
	}
	fprintf(fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			map.QNAME, map.FLAG, map.RNAME, map.POS+1, map.MAPQ, md_cg, map.MRNAME, map.MPOS, map.ISIZE, seq_cg, qual_cg);

	for ( i = 0; i < map.optSize; i++) {
		switch (map.optFields[i].type) {
			case 'A':
				fprintf(fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
				break;
			case 'i':
				fprintf(fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
				break;
			case 'f':
				fprintf(fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
				break;
			case 'Z':
			case 'H':
				fprintf(fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
				break;
		}
	}
	fprintf(fp, "\tGC:Z:%s", gc_cg);
	fprintf(fp, "\tGS:Z:%s", gs_cg);
	fprintf(fp, "\tGQ:Z:%s", gq_cg);
	if (readGroup[0] != 0)
		fprintf(fp, "\t%s:%c:%s", "RG", 'Z', readGroup);
	fprintf(fp, "\n");
}

void outputQ(SAM map) {
	int i = 0;
	int index = 0;
	int gap1 = 0;
	int gap2 = 0;
	int gap3 = 0;

	char gc_cg[SEQ_LENGTH];
	char md_cg[SEQ_LENGTH];
	char seq_cg[SEQ_LENGTH];
	char gs_cg[SEQ_LENGTH];
	char qual_cg[SEQ_LENGTH];
	char gq_cg[SEQ_LENGTH];

	if (readFormat == 1) {
		if(sscanf(map.CIGAR, "5M%dS10M%dN10M%dN10M", &gap1, &gap2, &gap3) == 3) {
			if(gap2 == 0)
				sprintf(md_cg, "%dM%dN10M", 25 - gap1, gap3);
			else
				sprintf(md_cg, "%dM%dN10M%dN10M", 15-gap1,gap2, gap3);
			sprintf(gc_cg, "%dS%dG%dS", 5-gap1, gap1, 30-gap1);
			for(i = 0; i < SEQ_LENGTH; i++) {
				if(i < 5 || i >= 5+gap1) {
					seq_cg[index] = map.SEQ[i];
					qual_cg[index] = map.QUAL[i];
					index++;
				}
			}
			seq_cg[index]='\0';
			qual_cg[index]='\0';
			index = 0;
			for(i = 0; i < SEQ_LENGTH; i++) {
				if(i >= 5-gap1 && i < 5+gap1) {
					gs_cg[index] = map.SEQ[i];
					gq_cg[index] = map.QUAL[i];
					index++;
				}
			}
			gs_cg[index]='\0';
			gq_cg[index]='\0';
		}
		else {
			sscanf(map.CIGAR, "10M%dN10M%dN10M%dS5M", &gap1, &gap2, &gap3);
			if(gap2 == 0)
				sprintf(md_cg, "10M%dN%dM", gap1, 25 - gap3);
			else
				sprintf(md_cg, "10M%dN10M%dN%dM", gap1, gap2, 15 - gap3);
			sprintf(gc_cg, "%dS%dG%dS", 30 - gap3, gap3, 5 - gap3);
			for(i = 0; i < SEQ_LENGTH; i++) {
				if(i < 30 || i >= 30 + gap3) {
					seq_cg[index] = map.SEQ[i];
					qual_cg[index] = map.QUAL[i];
					index++;
				}
			}
			seq_cg[index]='\0';
			qual_cg[index]='\0';
			index = 0;
			for(i = 0; i < SEQ_LENGTH; i++) {
				if(i >= 30 - gap3 && i < 30 + gap3) {
					gs_cg[index] = map.SEQ[i];
					gq_cg[index] = map.QUAL[i];
					index++;
				}
			}
			gs_cg[index]='\0';
			gq_cg[index]='\0';
		}
	}
	else {
		sscanf(map.CIGAR, "10M%dS10M%dS10M", &gap1, &gap2);
		sprintf(md_cg, "%dM", 30 - gap1 - gap2);
		sprintf(gc_cg, "%dS%dG%dS%dG%dS", 10, gap1, 10 - gap1, gap2, 10 - gap2);

		for(i = 0; i < SEQ_LENGTH; i++) {
			if((i < 10) ||
			 	 (i >= 10 + gap1 && i < 20 + gap1) ||
				 (i >= 20 + gap1 + gap2 && i < 30 + gap1 + gap2)) {
				seq_cg[index] = map.SEQ[i];
				qual_cg[index] = map.QUAL[i];
				index++;
			}
		}
		seq_cg[index]='\0';
		qual_cg[index]='\0';
		index = 0;
		for(i = 0; i < SEQ_LENGTH; i++) {
			if ((i >= 10 - gap1 && i < 10 + gap1) ||
					(i >= 20 - gap1 - gap2 && i < 20 + gap1 + gap2)){
				gs_cg[index] = map.SEQ[i];
				gq_cg[index] = map.QUAL[i];
				index++;
			}
		}
		gap3 = 0;
		gs_cg[index]='\0';
		gq_cg[index]='\0';
	}

	fprintf(_out_fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			map.QNAME, map.FLAG, map.RNAME, map.POS+1, map.MAPQ, md_cg, map.MRNAME, map.MPOS, map.ISIZE, seq_cg, qual_cg);
	for ( i = 0; i < map.optSize; i++) {
		switch (map.optFields[i].type) {
			case 'A':
				fprintf(_out_fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
				break;
			case 'i':
				fprintf(_out_fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
				break;
			case 'f':
				fprintf(_out_fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
				break;
			case 'Z':
			case 'H':
				fprintf(_out_fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
				break;
		}
	}
	fprintf(_out_fp, "\tGC:Z:%s", gc_cg);
	fprintf(_out_fp, "\tGS:Z:%s", gs_cg);
	fprintf(_out_fp, "\tGQ:Z:%s", gq_cg);
	if (readGroup[0] != 0)
		fprintf(_out_fp, "\t%s:%c:%s", "RG", 'Z', readGroup);
	fprintf(_out_fp, "\n");
}

int initOutput ( char *fileName, int compressed) {
  if (compressed) {
	  char newFileName[strlen(fileName)+4];
	  sprintf(newFileName, "%s.gz", fileName);
	  _out_gzfp = fileOpenGZ(newFileName, "w1f");
	  if (_out_gzfp == Z_NULL) {
	  	return 0;
		}
  	finalizeOutput = &finalizeGZOutput;
  	output = &gzOutputQ;
  	SAMheaderGZ(_out_gzfp);
	}
  else {
	  _out_fp = fileOpen(fileName, "w");
	  if (_out_fp == NULL) {
		  return 0;
		}
	  finalizeOutput = &finalizeTXOutput;
	  output = &outputQ;
	  SAMheaderTX(_out_fp, 1);
	}
  buffer[0] = '\0';
  return 1;
}

FILE * getOutputFILE() {
  if(_out_fp != NULL)
		return _out_fp;
//  else if(_out_gzfp != NULL)
//		return _out_gzfp;
  else
		return NULL;
}

void SAMheaderTX(FILE *outfp, int check) {
  FILE *fp;
  char fainame[FILE_NAME_LENGTH];
  char chrom[FILE_NAME_LENGTH];
  int chromlen;
  char rest[FILE_NAME_LENGTH];
  char *ret;

  sprintf(fainame, "%s.fai",fileName[0]);
  fp = fopen(fainame, "r");

  if (fp != NULL){
		fprintf(outfp, "@HD\tVN:1.4\tSO:unknown\n");

		while (fscanf(fp, "%s\t%d\t", chrom, &chromlen) > 0) {
		  ret = fgets(rest, FILE_NAME_LENGTH, fp);
		  fprintf(outfp, "@SQ\tSN:%s\tLN:%d\n", chrom, chromlen);
		}
		fclose(fp);

		if (readGroup[0] != 0 && sampleName[0] != 0)
			fprintf(outfp, "@RG\tID:%s\tSM:%s\tLB:%s\tPL:cg\n", readGroup, sampleName, libName);
		fprintf(outfp, "@PG\tID:sirFAST\tPN:sirFAST\tVN:%s.%s\n",  versionNumber, versionNumberF);
  }
  else if (check){
		fprintf(stdout, "WARNING: %s.fai not found, the SAM file(s) will not have a header.\n", fileName[0]);
		fprintf(stdout, "You can generate the .fai file using samtools. Please place it in the same directory with the index to enable SAM headers.\n");
  }

  ret++;
}

void SAMheaderGZ(gzFile outgzfp) {
  FILE *fp;
  char fainame[FILE_NAME_LENGTH];
  char chrom[FILE_NAME_LENGTH];
  int chromlen;
  char rest[FILE_NAME_LENGTH];
  char *ret;

  sprintf(fainame, "%s.fai",fileName[0]);
  fp = fopen(fainame, "r");

  if (fp != NULL){
	gzprintf(outgzfp, "@HD\tVN:1.4\tSO:unknown\n");

	while (fscanf(fp, "%s\t%d\t", chrom, &chromlen) > 0){
	  ret = fgets(rest, FILE_NAME_LENGTH, fp);
	  gzprintf(outgzfp, "@SQ\tSN:%s\tLN:%d\n", chrom, chromlen);
	}
	fclose(fp);

	if (readGroup[0] != 0 && sampleName[0] != 0)
	  gzprintf(outgzfp, "@RG\tID:%s\tSM:%s\tLB:%s\tPL:cg\n", readGroup, sampleName, libName);

	gzprintf(outgzfp, "@PG\tID:sirFAST\tPN:sirFAST\tVN:%s.%s\n",  versionNumber, versionNumberF);
  }
  else{
	fprintf(stdout, "WARNING: %s.fai not found, the SAM file(s) will not have a header.\n", fileName[0]);
	fprintf(stdout, "You can generate the .fai file using samtools. Please place it in the same directory with the index to enable SAM headers.\n");
  }

  ret++;
}
