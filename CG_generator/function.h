
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>

using namespace std;

int countRefSeq(string ref_filename, char**, int*, int*, int);
int readInput(string reads_file, string result_file, int margin);
void makeRefSeq(string ref_filename, char** ch_array, char** ch_array_error);
int rand_1st ();
int rand_2nd ();
int rand_3rd ();
int rand_4th ();
int rand_loc (int loc);
int addErr (char * ch_error, char ** ch_array, int ix, int jx, int no_error, int * storage);
int randErr (int width, int number, int * storage);
int randErrSeq ();
string getRefSeq(int coordinate, int size, string ref_filename);
struct expected_loc {
	int forward_loc;
	int reverse_loc;
	int forward_hit;
	int reverse_hit;
};
