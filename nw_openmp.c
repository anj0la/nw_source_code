/**
 * 
 * File: nw_openmp.c
 * Author: Anjola Aina
 * Last Modified: Thursday, April 18th, 2024
 * 
 * Information:
 * 
 * This function implements the Needleman-Wunsch algorithm in parallel using OpenMP directives. It uses a Cell matrix, holding the score and the direction of the sequence to make the traceback portion of the 
 * algorithm easier to compute.
 * 
 * For simplicity, we assume that the two sequences are of equal length.
 * 
 * Running the Program:
 * 
 * To run the program, specify the correct number of parameters (7), which are the following:
 *  - text1.txt - the first text file to write the first sequence into
 *  - text2.txt - the second text file to write the second sequence into
 *  - 4 - the length of both sequences
 *  - 8 - the number of threads
 *  - 1 - match score
 *  - -1 - mistmatch score
 *  - -2 - gap penalty
 * 
 * The code should be executed in your chosen terminal (system terminal or using an IDE for C: ./[executable_file_name] [first file] [second file] [length] [num_threads] [match] [mismatch] [gap_penalty]
 * Say we have our executable file compiled called nw_openmp, with two text files seq1.txt and seq2.txt and their respective lengths 4 characters, and a match, mismatch and gap penalty score of 1, -1 and -2 respectively.
 * Then, we would write the following line to the terminal:
 * 
 *      ./nw_openmp seq1.txt seq2.txt 4 8 1 -1 -2
 * 
 * Sources: 
 * 
 *  - strrev does not exist in the Linux implementation of C in <string.h>, so we used the following stack overflow resource to get the implementation of the funnction to save time: https://stackoverflow.com/questions/8534274/is-the-strrev-function-not-available-in-linux
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

// CONSTANTS
#define NUM_ARGS 7
#define TRUE 1
#define FALSE 0

// =============== STRUCTURES ===============

typedef enum Direction {
    DIAGONAL,
    UP,
    LEFT
} Direction;

typedef struct Cell {
    int score;
    Direction direction;
} Cell;

typedef struct Aligned_Sequences {
    char* aligned_s1;
    char* aligned_s2;
} Aligned_Sequences;

// =============== FUNCTION DEFINITIONS ===============

int max(int a, int b);
int min (int a, int b);
int check_if_correct_num_params_specified(int argc);
Cell max_cell(Cell a, Cell b, Cell c);
char* strrev(char* str);
void generate_rand_seqs(const char* file_1, const char* file_2, int seq_len);
char* get_seq(const char* filename);
Cell** allocate_matrix(int m, int n);
void init_matrix_p(Cell** matrix, int m, int n, int gap);
void fill_matrix_p(char* s1, char* s2, Cell** matrix, int thread_count, int match, int mis_match, int gap);
void traceback_matrix(char* s1, char* s2, Cell** matrix);
void nw_p(char* s1, char* s2, int thread_count, int match, int mis_match, int gap);

// =============== CODE ===============

int main(int argc, char* argv[]) {
    srand(time(NULL));
    // error checking for parameters
    int check_params = check_if_correct_num_params_specified(argc);
    if (check_params != TRUE) {
        return 1; // stop execution
    }

    generate_rand_seqs(argv[1], argv[2], atoi(argv[3]));

    char* s1 = get_seq(argv[1]);
	char* s2 = get_seq(argv[2]);

    // error checking for memory allocation
    if (s1 == NULL || s2 == NULL) {
        return 1; // stop execution
    }
    int thread_count = atoi(argv[4]);

    printf("thread count: %d\n", thread_count);
    
    int match = atoi(argv[5]);
    int mis_match = atoi(argv[6]);
    int gap = atoi(argv[7]);
    struct timeval start;
	struct timeval end;

    // timing program
    gettimeofday(&start, NULL);
    nw_p(s1, s2, thread_count, match, mis_match, gap); 
    gettimeofday(&end, NULL);

    // getting and printing result
    double seconds_s = (double) (end.tv_usec - start.tv_usec) / 1000000 +
         (double) (end.tv_sec - start.tv_sec);
    printf("Time taken to execute nw: %f seconds\n", seconds_s);

    free(s1);
    free(s2);

    return 0;
} // main

int max(int a, int b) {
    return a >= b ? a : b;
} // max

int min(int a, int b) {
    return a <= b ? a : b;
} // min

/**
 * Verifies that the correct number of parameters have been specified for the program.
 * @param argc number of arguments specified on the command line (excluding the executable file itself)
 * @return - 1 (true) if enough params have been specified, 0 (false) otherwise
 * 
*/
int check_if_correct_num_params_specified(int argc) {
    if (argc < NUM_ARGS || argc > NUM_ARGS) {
		return FALSE; // not enough or too much arguments were specified
	}
	return TRUE; // correct num of arguments were specified
} // check_if_correct_num_params_specified

/**
 * Returns the maximum value of the three specified values.
 * @param a the first value to be compared
 * @param b the second value to be compared
 * @param c the third value to be compared
 * @return - the maximum value from the three values
*/
Cell max_cell(Cell a, Cell b, Cell c) {
    Cell max_val = a;
    if (b.score > max_val.score) {
        max_val = b;
    }
    if (c.score > max_val.score) {
        max_val = c;
    }
    return max_val;
} // max

/**
 * Reverses the given string.
 * @param str the string to be reversed
 * @return - the string in its reversed form
*/
char* strrev(char *str) {
    if (!str || ! *str)
        return str;

    int i = strlen(str) - 1, j = 0;

    char ch;
    while (i > j) {
        ch = str[i];
        str[i] = str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
} // strrev

/**
Generates two random sequences of the specified length and places them in the specified text files.
@param file_1 the file that will contain the first sequence
@param file_2 the file that will contain the second sequence
@param seq_len - the length of each sequence
*/
void generate_rand_seqs(const char* file_1, const char* file_2, int seq_len) {
	FILE* fp_1 = fopen(file_1, "w");
	FILE* fp_2 = fopen(file_2, "w");
	char ch_dna[] = {'A', 'C', 'G', 'T'};
	if (fp_1 == NULL && fp_2 == NULL) {
		printf("Both files cannot be created and/or opened. \n");
		exit(1); // error occurred, stop execution of program
	} // if
	for (int i = 0; i < seq_len; i++) {
		int ch_seq_1_index = rand() % 4;
		int ch_seq_2_index = rand() % 4;
		fputc(ch_dna[ch_seq_1_index], fp_1);
		fputc(ch_dna[ch_seq_2_index], fp_2);
	} // for
	fclose(fp_1);
	fclose(fp_2);
} // generate_rand_seqs

/**
 * Gets and returns a sequence of characters from a file.
 * @param filename the file containing a sequence of characters
 * @return - a string containing the sequence of characters from the file
*/
char* get_seq(const char* filename) {
	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
        return NULL; // indicates that the file was null and no reading was achieved
    }

    fseek(fp, 0L, SEEK_END);
	int seq_len = ftell(fp);
    fseek(fp, 0L, SEEK_SET); // setting fp pointer back to the start when we read from the file again

	char* seq = (char*) malloc((seq_len + 1) * sizeof(char));
    fgets(seq, seq_len + 1, fp);
	fclose(fp);
	return seq;
} // get_seq

/**
 * Dynamically allocates a  reverses the string at the end matrix of size m and n, which are the lengths of the first and second sequence respectively.
 * @param m the length of the first sequence
 * @param n the length of the second sequence
 * @return - the cell matrix of size m x n
*/
Cell** allocate_matrix(int m, int n) {
    Cell** matrix = (Cell**) malloc((m + 1) * sizeof(Cell*));
    for (int i = 0; i < m + 1; i++) {
        matrix[i] = (Cell*) malloc((n + 1) * sizeof(Cell));
    }
    return matrix;
} // allocate_matrix

/**
 * Implements the initialization step of the algorithm, by initializing the first row and column of the cell matrix.
 * NOTE: This function is executed in parallel,
 * @param matrix the cell matrix to contain the all the scores and directions
 * @param size the length of the two sequences
 * @param gap the gap penalty score (determined on the command line)
*/
void init_matrix_p(Cell** matrix, int m, int n, int gap) {
    #pragma omp parallel for
    for (int i = 0; i < m + 1; i++) {
        matrix[i][0].score = i * gap;
    }
    #pragma omp parallel for
    for (int j = 0; j < n + 1; j++) {
        matrix[0][j].score = j * gap;
    }
} // init_matrix_p

/**
 * Implements the matrix filling step of the algorithm, by filling in all of the cells in the matrix using the anti-diagonal approach, 
 * to execute the algorithm in parallel, avoiding some of the data dependencies.
 * @param s1 the first sequence
 * @param s2 the second sequence
 * @param matrix the cell matrix to be filled wiih the scores and the matrix
 * @param match the match score (determined on the command line)
 * @param mis_match the mis_match score (determined on the command line)
 * @param gap the gap penalty (determined on the command line)
*/
void fill_matrix_p(char* s1, char* s2, Cell** matrix, int thread_count, int match, int mis_match, int gap) {
    int i, j, diag;
    int m = strlen(s1);
    int n = strlen(s2);
    int match_or_mismatch = 0;
    int gap_in_s1 = 0;
    int gap_in_s2 = 0;

    #pragma omp parallel num_threads(thread_count) default(none) shared(matrix, m, n, s1, s2, match, mis_match, gap) private(i, j, diag, match_or_mismatch, gap_in_s1, gap_in_s2) 
    {
    for (diag = 1; diag <= m + n - 1; diag++) {
        int start = max(1, diag - n + 1);
        int end = min(m, diag);
        #pragma omp for
        for (i = start; i <= end; i++) {
            j = diag - i + 1;
            match_or_mismatch = matrix[i - 1][j - 1].score + (s1[i - 1] == s2[j - 1] ? match : mis_match);
            gap_in_s2 = matrix[i - 1][j].score + gap;
            gap_in_s1 = matrix[i][j - 1].score + gap;

            Cell diagonal = {match_or_mismatch, DIAGONAL};
            Cell up = {gap_in_s2, UP};
            Cell left = {gap_in_s1, LEFT};

            Cell max_score = max_cell(diagonal, up, left);
            matrix[i][j].score = max_score.score;
            matrix[i][j].direction = max_score.direction;   
        } // inner for

        // synchronizes threads after completing computation for current anti-diagonal
        #pragma omp barrier
    } // outer for
} // parallel region
} // fill_matrix_p

/**
 * Implements the traceback step of the algorithm, by tracing back starting from the bottom right corner to get the optimal aligned sequences.
 * @param s1 the first sequence
 * @param s2 the second sequence
 * @param matrix the cell matrix (now filled)
*/
void traceback_matrix(char* s1, char* s2, Cell** matrix) {
    int i = strlen(s1);
    int j = strlen(s2);
    char* s1_align = (char*) malloc((i * j + 1) * sizeof(char));
    char* s2_align = (char*) malloc((j * j + 1) * sizeof(char));

    printf("Score: %d\n", matrix[i][j]);

    int k = 0;

    while (i > 0 || j > 0) {
        if (matrix[i][j].direction == DIAGONAL) {
            s1_align[k] = s1[i - 1];
            s2_align[k] = s2[j - 1];
            i--;
            j--;
        } else if (matrix[i][j].direction == UP) {
            s1_align[k] = s1[i - 1];
            s2_align[k] = '-';
            i--;
        } else {
            s1_align[k] = '-';
            s2_align[k] = s2[j - 1];
            j--;
        }
        k++;
    }
    s1_align[k] = '\0';
    s2_align[k] = '\0';

    // step 4 - reverse them and obtain the sequences 

    // print the aligned sequences (gets too long when sequence length is large)
    //printf("S1 Aligned: %s\n", strrev(s1_align));
    //printf("S2 Aligned: %s\n", strrev(s2_align));

    free(s1_align);
    free(s2_align);
} // traceback_matrix

/**
 * Implements the Needleman-Wunsch algorithm. It uses the following methods, which are
 * the required steps in the algorithm:
 *       init_matrix_p(matrix, m, n, gap);
 *       fill_matrix_p(s1, s2, matrix, match, mis_match, gap);
 *       traceback_matrix(s1, s2, matrix);
 * @param s1 the first sequence
 * @param s2 the second sequence
 * @param match the match score (determined on the command line)
 * @param mis_match the mis_match score (determined on the command line)
 * @param gap the gap penalty (determined on the command line)
*/
void nw_p(char* s1, char* s2, int thread_count, int match, int mis_match, int gap) {
    int m = strlen(s1);
    int n = strlen(s2);

    // allocate memory for cell matrix (bring into own function)
    Cell** matrix = allocate_matrix(m, n);

    // step 1- initalize matrix
    init_matrix_p(matrix, m, n, gap);

    // step 2 - fill in the matrix (DATA DEPENDECY)
    fill_matrix_p(s1, s2, matrix, thread_count, match, mis_match, gap);

    // step 3 - traceback matrix to get the aligned sequences (also has data dependecy i think)
    traceback_matrix(s1, s2, matrix);

    // freeing the matrix
    for (int i = 0; i < m + 1; i++) {
        free(matrix[i]);
    }
   free(matrix);
} // nw_p
