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
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include "Common.h"
#include "CommandLineParser.h"

int				uniqueMode=1;
int				indexingMode;
int				searchingMode;
int				pairedEndMode;
int				pairedEndModeMP;
int				pairedEndModePE;
int				pairedEndDiscordantMode;
int				transChromosomal=0;
int				pairedEndProfilingMode;
int				seqCompressed;
int				outCompressed;
int				cropSize = 0;
int				progressRep = 0;
int				minPairEndedDistance=-1;
int				maxPairEndedDistance=-1;
int				minPairEndedDiscordantDistance=-1;
int				maxPairEndedDiscordantDistance=-1;
int 				bestMode;
int 				nosamMode;
int                             debugMode=0;
char				*seqFile1;
char				*seqFile2;
char				*mappingOutput = "output";
char				*mappingOutputPath = "";
char				*unmappedOutput = "unmapped";
char				fileName[2][FILE_NAME_LENGTH];
int 				maxOEAOutput=100;
int 				maxDiscordantOutput=300;
unsigned char			errThreshold=255;
unsigned char			maxHits=0;
unsigned char			WINDOW_SIZE = 10;
unsigned int			CONTIG_SIZE;
unsigned int			CONTIG_MAX_SIZE;
char                            readGroup[FILE_NAME_LENGTH];
char                            sampleName[FILE_NAME_LENGTH];
char                            libName[FILE_NAME_LENGTH];
unsigned int			maxInputRead = -1;
int                             number_of_threads;

void printHelp();

int parseCommandLine (int argc, char *argv[])
{
  
  int o;
  int index;
  char *fastaFile = NULL;
  
  readGroup[0] = 0;
  sampleName[0] = 0;
  libName[0] = 0;
  number_of_threads = 1;

  static struct option longOptions[] = 
    {
      {"mp",   	        no_argument,	    &pairedEndModeMP,		1},
      {"pe",		no_argument,  	    &pairedEndModePE,		1},
      {"discordant-vh", no_argument,	    &pairedEndDiscordantMode,	1},
      {"trans",         no_argument,        &transChromosomal,          1},
      {"profile",       no_argument, 	    &pairedEndProfilingMode,	1},
      {"seqcomp",	no_argument,	    &seqCompressed,		1},
      {"outcomp",       no_argument,	    &outCompressed,		1},
      {"progress",	no_argument,	    &progressRep,		1},
      {"best",		no_argument,	    &bestMode,		1},
      {"debug",		no_argument,	    &debugMode,		1},
      {"index",		required_argument,  0, 			'i'},
      {"search",	required_argument,  0,			's'},
      {"help",		no_argument,	    0,			'h'},
      {"version",	no_argument,	    0,			'v'},
      {"seq",		required_argument,  0,			'x'},
      {"seq1",		required_argument,  0,			'x'},
      {"seq2",		required_argument,  0,			'y'},
      {"bucket",	required_argument,  0,			'b'},
      //      {"ws",		required_argument,  0,			'w'},
      {"min",		required_argument,  0,			'l'},
      {"max",		required_argument,  0,			'm'},
      {"crop",		required_argument,  0,			'c'},
      {"maxoea",        required_argument,  0,                  'a'},
      {"maxdis",        required_argument,  0,                  'd'},
      {"rg",            required_argument,  0,                  'g'},
      {"sample",        required_argument,  0,                  'p'},
      {"lib",           required_argument,  0,                  'r'},
      {"nosam",         no_argument,        &nosamMode,         1},
      {0,  0,  0, 0},
    };

  if (argc == 1){
    printHelp();
    return 0;
  }

  while ( (o = getopt_long ( argc, argv, "hvn:e:o:u:i:s:x:y:l:m:c:a:d:g:p:r:t:b:", longOptions, &index)) != -1 )
    {
      switch (o)
	{
	case 'b':					// DHL slice
	  maxInputRead = atoi(optarg);
	  break;
	case 'a':
	  maxOEAOutput = atoi(optarg);
	  if (maxOEAOutput == 0)
	    maxOEAOutput = 100000;
	  break;
	case 'd':
	  maxDiscordantOutput = atoi(optarg);
	  if (maxDiscordantOutput == 0)
	    maxDiscordantOutput = 100000;
	  break;
	case 'i':
	  indexingMode = 1;
	  fastaFile = optarg;
	  break;
	case 's':
	  searchingMode = 1;
	  fastaFile = optarg;
	  break;
	case 'c': 
	  cropSize = atoi(optarg);
	  break;
	  /*
	case 'w':
	  WINDOW_SIZE = atoi(optarg);
	  break; */
	case 'x':
	  seqFile1 = optarg;
	  break;
	case 't':
	  number_of_threads = atoi(optarg);
	  break;
	case 'y':
	  seqFile2 = optarg;
	  break;
	case 'u':
	  unmappedOutput = optarg;
	  break;
	case 'o':
	  mappingOutput = getMem(FILE_NAME_LENGTH, "mappingOutput @parseCommandLine()");
	  mappingOutputPath = getMem(FILE_NAME_LENGTH, "mappingOutput @parseCommandLine()");
	  stripPath (optarg, &mappingOutputPath, &mappingOutput);
	  break;
	case 'n':
	  maxHits = atoi(optarg);
	  break;
	case 'e':
	  errThreshold = atoi(optarg);
	  break;
	case 'l':
	  minPairEndedDistance = atoi(optarg);
	  break;
	case 'm':
	  maxPairEndedDistance = atoi(optarg);
	  break;					
	case 'g':
	  strcpy(readGroup, optarg);
	  break;					
	case 'p':
	  strcpy(sampleName, optarg);
	  break;					
	case 'r':
	  strcpy(libName, optarg);
	  break;					
	case 'h':
	  printHelp();
	  return 0;
	  break;
	case 'v':
	  fprintf(stdout, "sirFAST %s.%s with FastHASH\n", versionNumber, versionNumberF);
	  return 0;
	  break;				
	case '?': 
	  fprintf(stderr, "Unknown parameter: %s\n", longOptions[index].name);
	  abort();
	  return 0;
	  break;
	}

    }

  if (indexingMode + searchingMode != 1)
    {
      fprintf(stdout, "ERROR: Indexing / Searching mode should be selected\n");
      return 0;
    }

  //if (WINDOW_SIZE > 15 || WINDOW_SIZE < 11)
  /*
  if (WINDOW_SIZE > 15 || WINDOW_SIZE < 9)		// sirFAST
    {
      fprintf(stdout, "ERROR: Window size should be in [12..15]\n");
      return 0;
    }
  */

  if ( indexingMode )
    {
      CONTIG_SIZE	= 15000000;
      CONTIG_MAX_SIZE	= 40000000;

      if (fastaFile == NULL)
	{
	  fprintf(stdout, "ERROR: Reference(s) should be indicated for indexing\n");
	  return 0;
	}

      if (pairedEndDiscordantMode)
	{
	  fprintf(stdout, "ERROR: --discordant-vh cannot be used in indexing mode. \n");
	  return 0;
	}

    }


  if ( searchingMode )
    {
      CONTIG_SIZE	= 330000000;	// 2.5.0.4
      CONTIG_MAX_SIZE	= 330000000;	// 2.5.0.4
      
      if (pairedEndModeMP && pairedEndModePE){
	fprintf(stdout, "ERROR: Use either --pe or --mp, not both. \n");
	return 0;
      }
      
      pairedEndMode = (pairedEndModeMP || pairedEndModePE);

      if (maxInputRead == -1){
	if (pairedEndMode) maxInputRead = 250000;
	else maxInputRead = 5000000;
      }
  
		
      if (fastaFile == NULL)
	{
	  fprintf(stdout, "ERROR: Index File should be indiciated for searching\n");
	  return 0;
	}

      if (seqFile1 == NULL && seqFile2 == NULL)
	{
	  fprintf(stdout, "ERROR: Please indicate a sequence file for searching.\n");
	  return 0;
	}


      if (!pairedEndMode && seqFile2 != NULL)
	{
	  fprintf(stdout, "ERROR: Second File can be indicated in pairedend mode\n");
	  return 0;
	}

      if (pairedEndMode && (minPairEndedDistance <0 || maxPairEndedDistance < 0 || minPairEndedDistance > maxPairEndedDistance))
	{
	  fprintf(stdout, "ERROR: Please enter a valid range for pairedend sequences.\n");
	  return 0;
	}

      if (pairedEndMode && seqFile1 == NULL)
	{
	  fprintf(stdout, "ERROR: Please indicate the first file for pairedend search.\n");
	  return 0;
	}

      if (!pairedEndMode && pairedEndDiscordantMode)
	{
	  fprintf(stdout, "ERROR: --discordant-vh should be used with --pe\n");
	  return 0;
	}

      if (!pairedEndMode && pairedEndProfilingMode)
	{
	  fprintf(stdout, "ERROR: --profile should be used with --pe\n");
	  return 0;
	}

      if (pairedEndMode)
	pairedEndDiscordantMode = 1;
		
      if (readGroup[0] != 0 && sampleName[0] == 0)
	{
	  fprintf(stdout, "ERROR: --sample should be used with --rg\n");
	  return 0;
	}

      if (readGroup[0] != 0 && libName[0] == 0)
	{
	  fprintf(stdout, "ERROR: --lib should be used with --rg\n");
	  return 0;
	}

      if (readGroup[0] == 0 && sampleName[0] != 0)
	{
	  fprintf(stdout, "ERROR: --rg should be used with --sample\n");
	  return 0;
	}

      if (readGroup[0] == 0 && libName[0] != 0)
	{
	  fprintf(stdout, "ERROR: --rg should be used with --lib\n");
	  return 0;
	}

    }

  sprintf(fileName[0], "%s", fastaFile);
  sprintf(fileName[1], "%s.index", fileName[0]); 
       


  if (pairedEndProfilingMode)
    {

      minPairEndedDistance = 0;
      maxPairEndedDistance = 300000000;

    }

  if (pairedEndDiscordantMode)
    {
      minPairEndedDiscordantDistance = minPairEndedDistance;
      maxPairEndedDiscordantDistance = maxPairEndedDistance;

      minPairEndedDistance = 0;
      maxPairEndedDistance = 300000000;
    }
  
  return 1;
}


void printHelp()
{
  char *errorType;

  fprintf(stdout,"sirFAST : Short interrupted Read Fast Alignment Search Tool. Enhanced with FastHASH.\n\n");
  fprintf(stdout,"Usage: sirfast [options]\n\n");
  errorType="hamming distance";
  

  fprintf(stdout,"General Options:\n");
  fprintf(stdout," -v|--version\t\tCurrent Version.\n");
  fprintf(stdout," -h\t\t\tShows the help file.\n");
  fprintf(stdout,"\n\n");

  fprintf(stdout,"Indexing Options:\n");
  fprintf(stdout," --index [file]\t\tGenerate an index from the specified fasta file. \n");
  //fprintf(stdout," --ws [int]\t\tSet window size for indexing (default:12 max:14).\n");
  fprintf(stdout,"\n\n");

  fprintf(stdout,"Searching Options:\n");
  fprintf(stdout," --search [file]\tSearch in the specified genome. Provide the path to the fasta file. \n\t\t\tIndex file should be in the same directory.\n");
  fprintf(stdout," --pe \t\t\tSearch will be done in Paired-End mode.\n");
  fprintf(stdout," --seq [file]\t\tInput sequences in fasta/fastq format [file]. If \n\t\t\tpaired end reads are interleaved, use this option.\n");
  fprintf(stdout," --seq1 [file]\t\tInput sequences in fasta/fastq format [file] (First \n\t\t\tfile). Use this option to indicate the first file of \n\t\t\tpaired end reads. \n");
  fprintf(stdout," --seq2 [file]\t\tInput sequences in fasta/fastq format [file] (Second \n\t\t\tfile). Use this option to indicate the second file of \n\t\t\tpaired end reads.  \n");
  fprintf(stdout," -o [file]\t\tOutput of the mapped sequences. The default is \"output\".\n");
  fprintf(stdout," -u [file]\t\tSave unmapped sequences in fasta/fastq format.\n");
  fprintf(stdout," --best   \t\tOnly the best mapping from all the possible mapping is returned.\n");
  fprintf(stdout," --seqcomp \t\tIndicates that the input sequences are compressed (gz).\n");
  fprintf(stdout," --outcomp \t\tIndicates that output file should be compressed (gz).\n");
  fprintf(stdout," -e [int]\t\tMaximum allowed %s (default 4%% of the read length).\n", errorType);
  fprintf(stdout," --min [int]\t\tMin distance allowed between a pair of end sequences.\n");
  fprintf(stdout," --max [int]\t\tMax distance allowed between a pair of end sequences.\n");
  fprintf(stdout," --maxoea [int]\t\tMax number of One End Anchored (OEA) returned for each read pair.\n\t\t\tWe recommend 100 or above for NovelSeq use. Default = 100.\n");
  fprintf(stdout," --maxdis [int]\t\tMax number of discordant map locations returned for each read pair.\n\t\t\tWe recommend 300 or above for VariationHunter use. Default = 300.\n");
  fprintf(stdout," --crop [int]\t\tTrim the reads to the given length.\n");
  fprintf(stdout," --sample [string]\tSample name to be added to the SAM header (optional).\n");
  fprintf(stdout," --rg [string]\t\tRead group ID to be added to the SAM header (optional).\n");
  fprintf(stdout," --lib [string]\t\tLibrary name to be added to the SAM header (optional).\n");
  fprintf(stdout,"\n\n");
}
