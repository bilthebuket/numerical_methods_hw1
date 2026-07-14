#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct Matrix
{
	double* vals;
	int rows;
	int cols;
} Matrix;

// creates a matrix from a string
// format (first number is row, second number is column): 11, 12, 13; 21, 22, 23; 31, 32, 33
Matrix* matrix_create(char* s);
Matrix* matrix_duplicate(Matrix* m);
void matrix_free(Matrix* m);

// modifies add_to and sub_from
void matrix_add_matrix(Matrix* add_to, Matrix* add_from);
void matrix_sub_matrix(Matrix* sub_from, Matrix* to_sub);

// returns the multiplied matrix
Matrix* matrix_mult_matrix(Matrix* left, Matrix* right);
Matrix* matrix_mult_matrix_multithread(Matrix* left, Matrix* right);

// because i wanted to compare my optimizations 
Matrix* matrix_mult_matrix_bad(Matrix* left, Matrix* right);
Matrix* matrix_mult_matrix_decent(Matrix* left, Matrix* right);

// adds/subtracts/multiplies every value of the matrix by the constant
void matrix_add_constant(Matrix* m, double d);
void matrix_sub_constant(Matrix* m, double d);
void matrix_mult_constant(Matrix* m, double d);

void matrix_swap_row(Matrix* m, int row1, int row2);
// swaps is an integer array of the out of order rows
// ex: if swaps = 1, 0, 2, then in matrix m, the first and second rows will be swapped
void matrix_swap_rows(Matrix* m, int* swaps);

// for each column i, apply_to[i] += to_apply[i] * multiplier
void matrix_add_row(Matrix* m, int to_apply, int apply_to, double multiplier, bool truncate);
// multiplies every value in the row by val
void matrix_mult_row(Matrix* m, int row, double val, bool truncate);

// assumes m is a diagnol matrix that just needs to be reduced to the identity matrix by dividing each row by the value on the diagnol
void matrix_reduce_diagnol(Matrix* m, bool truncate);
// assumes below the diagnol is already zeros, is essentially back substitution
void matrix_reduce_upper(Matrix* m, bool truncate);

// both of these return int arrays which have the new ordering of the rows
int* matrix_reduce_lower(Matrix* m, bool partial_pivot, bool truncate);
// reduces m but leaves L of the LU decomposition below the diagnol instead of zeros
int* matrix_decompose(Matrix* m, bool partial_pivot, bool truncate);

// m becomes L where L * L_transpose = m
void matrix_cholesky_decompose(Matrix* m);

void matrix_invert(Matrix* m);
void matrix_transpose(Matrix* m);

// returns a new matrix with everything from m below the diagnol, zeros otherwise
Matrix* matrix_get_lower(Matrix* m);
// returns a new matrix with everything from m above the diagnol (does not include the diagnol itself), zeros otherwise
Matrix* matrix_get_upper(Matrix* m);
// returns a new matrix with everything from m on the diagnol, zeros otherwise
Matrix* matrix_get_diagnol(Matrix* m);

Matrix* matrix_make_identity(int size);

// augments add_from onto add_to, frees add_from
void matrix_augment(Matrix* add_to, Matrix* add_from);
// returned matrix contains columns including cutoff_column to the end, m contains the remaining columns
Matrix* matrix_unaugment(Matrix* m, int cutoff_column);

// sums sqaure of every value in m, then returns the square root of that sum
double matrix_get_magnitude(Matrix* m);

void matrix_print(Matrix* m, FILE* f);

bool matrix_equals(Matrix* m1, Matrix* m2);

#endif
