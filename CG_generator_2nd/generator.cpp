#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <random>
#include <iostream>
#include "function.h"

int main(int argc, char* argv[]) {
	if (argc < 4) {
		fprintf(stdout, "Need parameter [reference/error/size/fraction_of_error_read]\n");
		return 0;
	}

	int no_error = atoi(argv[2]);
	int total_reads = atoi(argv[3]);
	int fraction_of_error_read = atoi(argv[4]);
	int tmp_no_error = 0;

//---------- System Define ----------
//	int total_reads = 100000;
//	int no_error = 0;
//-----------------------------------

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(400, 10);

	char * ch_name_array[100];
	int ch_count_array[100];
	int reads_count[100];
	char * ch_array[100];
	char * ch_array_error[100];
	char ch_error[20];
	int gap_1, gap_2, gap_4, gap_6, gap_7;
	int base;
	int loc;
	int ErrSeqForward = 100;
	int ErrSeqReverse = 100;
	int storageForward[10];
	int storageReverse[10];

	for (int x = 0; x < 100; x++)
		ch_name_array[x] = (char*) malloc(sizeof(char)*100);
	int no_ch = countRefSeq(argv[1], ch_name_array, ch_count_array,
							reads_count, total_reads);
	for (int y = 0; y < no_ch; y++) {
		ch_array[y] = (char*) malloc(sizeof(char)*(ch_count_array[y]+1));
		ch_array_error[y] = (char*) malloc(sizeof(char)*(ch_count_array[y]+1));
	}

	makeRefSeq(argv[1], ch_array, ch_array_error);

	printf("#ASSEMBLY_ID    NA12878-200-37-ASM\n");
	printf("#BATCH_FILE_NUMBER  1\n");
	printf("#BATCH_OFFSET   0\n");
	printf("#FIELD_SIZE 460800\n");
	printf("#FORMAT_VERSION 2.0\n");
	printf("#GENERATED_AT   2011-Oct-08 00:33:17.041145\n");
	printf("#GENERATED_BY   exportTools\n");
	printf("#LANE   L01\n");
	printf("#LIBRARY    GS00664-CLS_G03\n");
	printf("#SAMPLE GS00362-DNA_C04\n");
	printf("#SLIDE  GS21184-FS3\n");
	printf("#SOFTWARE_VERSION   2.0.0.26\n");
	printf("#TYPE   READS\n");
	printf("\n");
	printf(">flags  reads   scores\n");

//---------- Local System Define ----------
//	no_ch = 1;
//	reads_count[0] = 10000;
//-----------------------------------------
	no_ch = 1;
	reads_count[0] = 1000;

	for (int j = 0; j < no_ch; j++) {
		for (int k = 0; k < reads_count[j]; k++) {
			gap_1 = rand_1st();
			gap_2 = rand_2nd();
			gap_4 = int(distribution(generator));
			gap_6 = rand_2nd();
			gap_7 = rand_1st();

			if (rand()%100 < fraction_of_error_read) {
				tmp_no_error = no_error;
			}
			else {
				tmp_no_error = 0;
			}

			if (tmp_no_error != 0) {
				ErrSeqForward = randErrSeq();
				ErrSeqReverse = randErrSeq();

//---------- Error Segment Definition ----------
//				ErrSeqForward = 3;
//				ErrSeqReverse = 1;
//----------------------------------------------
				if (ErrSeqForward == 0)
					randErr (10, tmp_no_error, storageForward);
				else
					randErr (10, tmp_no_error, storageForward);
				if (ErrSeqReverse == 0)
					randErr (10, tmp_no_error, storageReverse);
				else
					randErr (10, tmp_no_error, storageReverse);
			}

			base = rand_loc(ch_count_array[j]);
			if (base - 500 > ch_count_array[j])
				base = base - 500;
			loc = base;

			int detect_n = 0;
			for (int k = 0; k < 400; k++){
				if (ch_array[j][base + k] == 'N' || ch_array[j][base + k] == 'n')
					detect_n = 1;
			}

			if (detect_n == 1) {
				k--;
				continue;
			}

			for (int x = 0; x < 10; x++) {
				ch_error[x] = 'a';
			}

			printf("5	");

			//------------------------------------------------------------
			if (ErrSeqForward == 0) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageForward);
				for(int i = 0; i < 10; i++)
					printf("%c", ch_error[i]);
			}
			else {
				for(int i = loc; i < loc+10; i++)
					printf("%c", ch_array[j][i]);
			}
			//------------------------------------------------------------
			loc = loc + 10 + gap_1;
			if (ErrSeqForward == 1) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageForward);
				for(int i = 0; i < 9; i++)
					printf("%c", ch_error[i]);
				printf("N");
			}
			else {
				for(int i = loc; i < loc+9; i++)
					printf("%c", ch_array[j][i]);
				printf("N");

			}
			//------------------------------------------------------------
			loc = loc + 10 + gap_2;
			if (ErrSeqForward == 2) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageForward);
				for(int i = 0; i < 10; i++)
					printf("%c", ch_error[i]);
			}
			else {
				for(int i = loc; i < loc+10; i++)
					printf("%c", ch_array[j][i]);
			}
			//------------------------------------------------------------
			loc = loc + 10 + gap_4;
			if (ErrSeqReverse == 2) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageReverse);
				for(int i = 0; i < 10; i++)
					printf("%c", ch_error[i]);
			}
			else {
				for(int i = loc; i < loc+10; i++)
					printf("%c", ch_array[j][i]);
			}
			//------------------------------------------------------------
			loc = loc + 10 + gap_6;
			if (ErrSeqReverse == 1) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageReverse);
				printf("N");
				for(int i = 1; i < 10; i++)
					printf("%c", ch_error[i]);
			}
			else {
				printf("N");
				for(int i = loc + 1; i < loc+10; i++)
					printf("%c", ch_array[j][i]);
			}
			//------------------------------------------------------------
			loc = loc + 10 + gap_7;
			if (ErrSeqReverse == 0) {
				addErr (ch_error, ch_array, j, loc, tmp_no_error, storageReverse);
				for(int i = 0; i < 10; i++)
					printf("%c", ch_error[i]);
			}
			else {
				for(int i = loc; i < loc+10; i++)
					printf("%c", ch_array[j][i]);
			}
			//------------------------------------------------------------
			printf("	######################################################################");
			printf("	Loc: %d Gap: %d %d %d %d %d ",
					base, gap_1, gap_2, gap_4, gap_6, gap_7);
			printf(" Error: %d", tmp_no_error);
			printf(" Ch_Name: %s", ch_name_array[j]);
			//printf(" Ferr: %d %d", ErrSeqForward, storageForward[0]);
			//printf(" Rerr: %d %d", ErrSeqReverse, storageReverse[0]);
			//printf(" No_Err: %d", no_error);
		}
	}
	return 0;
}
