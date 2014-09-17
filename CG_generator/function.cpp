#include <iostream>
#include <fstream>
#include "function.h"

int addErr (char * ch_error, char ** ch_array, int ix, int jx, int no_error, int * storage) {

	for (int k = 0; k < 10; k++) {
		ch_error[k] = ch_array[ix][jx + k];
	}
	for (int i = 0; i < no_error; i++) {
		char tmp1 = ch_error[storage[i]];
		char tmp2;
		if (tmp1 == 'A' || tmp1 == 'a')
			//tmp2 = '1';
			tmp2 = 'C';
		else if (tmp1 == 'C' || tmp1 == 'c')
			//tmp2 = '1';
			tmp2 = 'G';
		else if (tmp1 == 'G' || tmp1 == 'g')
			//tmp2 = '1';
			tmp2 = 'T';
		else if (tmp1 == 'T' || tmp1 == 't')
			//tmp2 = '1';
			tmp2 = 'A';
		ch_error[storage[i]] = tmp2;
	}
	return 0;
}

int readInput(string reads_file, string result_file, int margin) {
	FILE *pFileRead;
	int loc, gap1, gap2, gap3, gap4, gap5, gap6, gap7;
	char s0[10], s1[100], s2[100], s3[10], s4[10], s5[10], s6[10];
	char tmp;
	int total_reads = 0;
	expected_loc * reads_loc;
	int order = 0;

	pFileRead  = fopen (reads_file.c_str(), "r");

	while (!feof(pFileRead)) {
		tmp = fgetc (pFileRead);
		if (feof(pFileRead))
			break;
		else if (tmp == '>' || tmp  == '#' || tmp == '\n')
			while (tmp != '\n')
				tmp = fgetc (pFileRead);
		else {
			fseek (pFileRead, -1, SEEK_CUR);
			fscanf (pFileRead, "%s %s %s %s %d %s %d %d %d %d %d %d %d %s %s",
				s0, s1, s2, s3, &loc, s4, &gap1, &gap2, &gap3, &gap4, &gap5, &gap6, &gap7, s5, s6);
			if (feof(pFileRead))
				break;
			else
				total_reads++;
		}
	}

	fseek (pFileRead, 0, SEEK_SET);
	reads_loc = (expected_loc*) malloc(sizeof(expected_loc)*total_reads);

	while (!feof(pFileRead)) {
		tmp = fgetc (pFileRead);
		if (feof(pFileRead))
			break;
		else if (tmp == '>' || tmp  == '#' || tmp == '\n')
			while (tmp != '\n')
				tmp = fgetc (pFileRead);
		else {
			fseek (pFileRead, -1, SEEK_CUR);
			fscanf (pFileRead, "%s %s %s %s %d %s %d %d %d %d %d %d %d %s %s",
					s0, s1, s2, s3, &loc, s4, &gap1, &gap2, &gap3, &gap4, &gap5, &gap6, &gap7, s5, s6);
			if (feof(pFileRead))
				break;
			reads_loc[order].forward_loc = loc;
			reads_loc[order].reverse_loc = loc + gap1 + gap2 + gap3 + gap4 + 1 + 35;
			reads_loc[order].forward_hit = 0;
			reads_loc[order].reverse_hit = 0;
			order++;
		}
	}
	fclose (pFileRead);

	FILE *pFileResult;
	pFileResult  = fopen (result_file.c_str(), "r");
	char t0[40],   t1[40], t2[40], t3[40], t4[40], t5[40], t6[40], t7[40], t8[40], t9[40];
	char t10[40], t11[40], t12[40], t13[40], t14[40], t15[40], t16[40], t17[40], t18[18];
	int  n1, n2, n3, n4;
	int out = 1;

	while (!feof(pFileResult)) {
		fscanf (pFileResult,
			"%[a-z1]%[_]%[a-z1]%[_]%[a-z1]%[_]%d%[/]%d %d %s %d %s %s %s %s %s %s %s %s %s \n",
			t0, t1, t2, t3, t4, t5, &n1, t6, &n2, &n3, t8, &n4,
			t9, t10, t11, t12, t13, t14, t15, t16, t17);
		if (n2 == 0) {
			if (reads_loc[n1].forward_loc >= n4 - margin && reads_loc[n1].forward_loc <= n4 + margin)
				reads_loc[n1].forward_hit = 1;
		}
		else {
			if (reads_loc[n1].reverse_loc >= n4 - margin && reads_loc[n1].reverse_loc <= n4 + margin)
				reads_loc[n1].reverse_hit = 1;
		}
	}
	fclose (pFileResult);

	int map_success_forward = 0;
	int map_success_reverse = 0;
	for (int i = 0; i < total_reads; i++) {

		cout << "##" << i << "## " << reads_loc[i].forward_loc << " " << reads_loc[i].forward_hit;
		cout << " -- " << reads_loc[i].reverse_loc << " " << reads_loc[i].reverse_hit << endl;

		if (reads_loc[i].forward_hit == 1)
			map_success_forward++;
		if (reads_loc[i].reverse_hit == 1)
			map_success_reverse++;
	}
	cout << endl << "Total Margin: " << margin;
	cout << endl << "Total Forward Hit: " << map_success_forward;
	cout << endl << "Total Reverse Hit: " << map_success_reverse << endl;
	return 0;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void makeRefSeq(string ref_file, char ** ch_array, char ** ch_array_error){
	FILE *pFile;
	int m = 0;
	int i = 0;
	int init = 0;
	int no_ch = 0;
	char temp;

	pFile = fopen (ref_file.c_str(), "r");

	while (!feof(pFile)) {
		temp = fgetc (pFile);
		if (temp == '>') {
			if (init == 0) {
				init = 1;
			}
			else {
				ch_array[no_ch][m] = '\0';
				ch_array_error[no_ch][m] = '\0';
				no_ch++;
			}
			while (temp != '\n'){
				temp = fgetc (pFile);
			}
			m = 0;
		}
		else if (temp != '\n') {
			ch_array[no_ch][m] = temp;
			ch_array_error[no_ch][m] = temp;
			m++;
		}
	}
	ch_array[no_ch][m-1] = '\0';
	ch_array_error[no_ch][m-1] = '\0';
	fclose (pFile);
}

int countRefSeq(string ref_file, char ** ch_name_array, int * ch_count_array,
				int * reads_count, int total_reads) {
	FILE *pFile;
	int n = 0;
	int m = 0;
	int i = 0;
	int init = 0;
	int no_ch = 0;
	char temp;
	char * ch_name;
	unsigned int total_ch = 0;

	ch_name = (char*) malloc(sizeof(char)*100);
	pFile = fopen (ref_file.c_str(), "r");

	while (!feof(pFile)) {
		temp = fgetc (pFile);
		n++;
		if (temp == '>') {

			if (init == 0)
				init = 1;
			else
				ch_count_array[no_ch-1] = m;

			i = 0;
			while (temp != '\n'){
				temp = fgetc (pFile);
				ch_name[i] = temp;
				ch_name_array[no_ch][i] = temp;
				i++;
			}
			ch_name[i] = '\0';
			ch_name_array[no_ch][i] = '\0';

			m = 0;
			no_ch++;
		}
		else if (temp != '\n')
			m++;
	}
	ch_count_array[no_ch-1] = m-1;
	fclose (pFile);

	i = 0;
	for (i = 0; i < no_ch; i++) {
		total_ch = total_ch + ch_count_array[i];
	}
	unsigned int test = 0;
	for (i = 0; i < no_ch; i++) {
		reads_count[i] = (int)(((double)ch_count_array[i]/(double)total_ch)*(double)total_reads);
		test = test + reads_count[i];
	}

	for (i = 0; i < total_reads - test; i++) {
		reads_count[i]++;
	}

	return no_ch;
}

int randErrSeq () {
	return rand() % 4;
}

int randErr (int width, int number, int * storage){
	int tmp = 0;
	for (int i = 0; i < number; ) {
		tmp = rand() % width;
		int flag_same = 0;
		for (int j = 0; j < i; j++) {
			if (storage[j] == tmp) {
				flag_same = 1;
			}
		}
		if (flag_same == 0) {
			storage[i] = tmp;
//			printf( " %c%d ", storage[i]);
			i++;
		}
	}
	// DEBUG------------------------------------------
//	printf("Storage: ");
//	for (int k = 0; k < number; k++) {
//		printf("%d ", storage[k]);
//	}
//	printf("\n");
	// DEBUG------------------------------------------
	return 0;
}

int rand_1st (){
	int v1 = rand() % 100;
	if (v1 <= 7) return -3;
	else if (v1 <= 14) return -1;
	else return -2;
}

int rand_2nd (){
	int v1 = rand() % 100;
	if (v1 <= 90) return 0;
	else if (v1 <= 97) return 1;
	else return 2;
}

int rand_3rd (){
	int v1 = rand() % 100;
	if (v1 <= 30) return 5;
	else if (v1 <= 90) return 6;
	else return 7;
}

int rand_loc (int loc){
	return rand() % loc;
}


string getRefSeq(int coordinate, int size,  string ref_filename) {
	int REF_TABLE_SIZE = 20;
    string result_string(100, 'A');
    FILE * pFileS;
    char * search_string;
    int  size_opt;
    int  coordinate_opt;
    int  boundary;
    int  boundary_detect = 0;
    int  read_number;

    //Initialize string size
    //result_string.resize(REF_TABLE_SIZE);

    boundary    = coordinate-(coordinate/REF_TABLE_SIZE)*REF_TABLE_SIZE + size; // boundary checking
    size_opt    = size + boundary / REF_TABLE_SIZE + 2;             // size optimization
    coordinate_opt  = coordinate+coordinate/REF_TABLE_SIZE;             // add # of new line characters
    search_string   = (char*) malloc(sizeof(char)*size_opt);
    pFileS = fopen (ref_filename.c_str(), "r");
    fseek (pFileS, 0, SEEK_SET );
    read_number = fread (search_string, 1, size_opt, pFileS);
    if (search_string[0] == '\n') {
        fseek (pFileS, coordinate_opt+1, SEEK_SET );
        read_number = fread (search_string, 1, size_opt, pFileS);
    } else {
        fseek (pFileS, coordinate_opt, SEEK_SET );
        read_number = fread (search_string, 1, size_opt, pFileS);
    }
    for (int i = 0; i < size ; i++) {
        if ((search_string[i+boundary_detect] == '\n')) {
            boundary_detect = boundary_detect + 1;
        }
        result_string[i] = search_string[i+boundary_detect];
        if (boundary_detect == 2){
//cout << " DETECTED 2 BOUNDARY " << boundary_detect << endl;
        }
    }
    fclose (pFileS);
    free(search_string);
    return result_string;
}
