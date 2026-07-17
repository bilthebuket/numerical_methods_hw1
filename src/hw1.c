#include <stdio.h>
#include "matrix.h"

// solves a system using gaussian elimination and lu decomposition
Matrix* solve_system(Matrix* m, bool partial_pivot, bool truncate)
{
	Matrix* b = matrix_unaugment(m, m->cols - 1);
	int* swaps = matrix_decompose(m, partial_pivot, truncate);
	matrix_swap_rows(b, swaps);
	Matrix* lower = matrix_get_lower(m);
	Matrix* upper = matrix_get_upper(m);
	Matrix* i = matrix_make_identity(lower->rows);
	Matrix* d = matrix_get_diagnol(m);
	matrix_add_matrix(lower, i);
	matrix_add_matrix(upper, d);
	matrix_free(m);
	matrix_free(i);
	matrix_free(d);

	Matrix* lower_dupe = matrix_duplicate(lower);
	Matrix* upper_dupe = matrix_duplicate(upper);

	matrix_augment(lower, b);
	int* arr = matrix_reduce_lower(lower, partial_pivot, truncate);
	if (arr != NULL)
	{
		free(arr);
	}
	matrix_reduce_diagnol(lower, truncate);
	b = matrix_unaugment(lower, lower->cols - 1);
	matrix_augment(upper, b);
	matrix_reduce_upper(upper, truncate, false);
	matrix_reduce_diagnol(upper, truncate);
	b = matrix_unaugment(upper, upper->cols - 1);

	printf("Solution:\n");
	matrix_print(b, stdout);
	printf("L:\n");
	matrix_print(lower_dupe, stdout);
	printf("U:\n");
	matrix_print(upper_dupe, stdout);

	Matrix* mult = matrix_mult_matrix(lower_dupe, upper_dupe);
	printf("LU:\n");
	matrix_print(mult, stdout);
	matrix_free(mult);

	matrix_free(upper);
	matrix_free(lower);
	matrix_free(lower_dupe);
	matrix_free(upper_dupe);
	free(swaps);
	return b;
}

static int problem_two(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("please provide exactly one system to solve and its solution\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}

	Matrix* computed = solve_system(m, false, false);
	Matrix* actual = matrix_create(argv[3]);
	Matrix* A = matrix_create(argv[2]);
	Matrix* b = matrix_unaugment(A, A->cols - 1);

	Matrix* mult = matrix_mult_matrix(A, computed);
	matrix_sub_matrix(b, mult);
	double error2 = matrix_get_magnitude(b);

	matrix_free(mult);
	matrix_free(b);
	matrix_free(A);

	matrix_sub_matrix(computed, actual);
	double error1 = matrix_get_magnitude(computed);

	matrix_free(computed);
	matrix_free(actual);

	printf("||x_computed - x_true||_2: %lf\n", error1);
	printf("||b - Ax_computed||_2: %lf\n", error2);

	return 0;
}

static int problem_three(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("please provide exactly one matrix to decompose\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}

	matrix_cholesky_decompose(m);
	Matrix* lower = matrix_get_lower(m);
	Matrix* d = matrix_get_diagnol(m);
	matrix_add_matrix(lower, d);
	matrix_free(d);
	matrix_free(m);
	Matrix* transpose = matrix_duplicate(lower);
	matrix_transpose(transpose);

	printf("L:\n");
	matrix_print(lower, stdout);
	printf("\nL_transpose:\n");
	matrix_print(transpose, stdout);
	matrix_free(lower);
	matrix_free(transpose);
	return 0;
}

static int problem_four(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("please provide exactly one system to solve and its solution\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}
	Matrix* actual = matrix_create(argv[3]);

	Matrix* no_swap = matrix_duplicate(m);
	int* arr = matrix_reduce_lower(no_swap, false, true);
	if (arr != NULL)
	{
		free(arr);
	}
	matrix_reduce_upper(no_swap, true, true);
	Matrix* b = matrix_unaugment(no_swap, no_swap->cols - 1);
	printf("truncated, no pivot:\n");
	matrix_print(b, stdout);

	Matrix* actual_dupe = matrix_duplicate(actual);
	matrix_sub_matrix(actual_dupe, b);
	double error = matrix_get_magnitude(actual_dupe) / matrix_get_magnitude(actual);
	printf("Relative Error: %lf\n", error);

	matrix_free(actual_dupe);
	matrix_free(b);
	matrix_free(no_swap);

	Matrix* swap = matrix_duplicate(m);
	arr = matrix_reduce_lower(swap, true, true);
	matrix_reduce_upper(swap, true, true);
	b = matrix_unaugment(swap, swap->cols - 1);
	if (arr != NULL)
	{
		free(arr);
	}
	printf("truncated, pivoted:\n");
	matrix_print(b, stdout);

	actual_dupe = matrix_duplicate(actual);
	matrix_sub_matrix(actual_dupe, b);
	error = matrix_get_magnitude(actual_dupe) / matrix_get_magnitude(actual);
	printf("Relative Error: %lf\n", error);

	matrix_free(actual_dupe);
	matrix_free(b);
	matrix_free(swap);

	printf("not truncated, pivoted, lu decomposition:\n");
	Matrix* A = matrix_duplicate(m);
	b = matrix_unaugment(A, A->cols - 1);
	Matrix* sol = solve_system(m, true, false);
	Matrix* mult = matrix_mult_matrix(A, sol);
	matrix_sub_matrix(b, mult);
	error = matrix_get_magnitude(b);
	printf("Residual Norm: %lf\n", error);
	matrix_free(sol);
	matrix_free(mult);
	matrix_free(A);
	matrix_free(b);
	matrix_free(actual);
	
	return 0;
}

static int problem_five(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("please provide the input data as one matrix\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}

	Matrix* b = matrix_unaugment(m, m->cols - 1);
	Matrix* transpose = matrix_duplicate(m);
	matrix_transpose(transpose);
	Matrix* M = matrix_mult_matrix(transpose, m);
	Matrix* c = matrix_mult_matrix(transpose, b);
	matrix_free(transpose);

	matrix_cholesky_decompose(M);
	Matrix* lower = matrix_get_lower(M);
	Matrix* d = matrix_get_diagnol(M);
	matrix_add_matrix(lower, d);
	matrix_free(d);
	matrix_free(M);
	transpose = matrix_duplicate(lower);
	matrix_transpose(transpose);

	matrix_augment(lower, c);
	int* arr = matrix_reduce_lower(lower, true, false);
	matrix_reduce_upper(lower, false, false);
	matrix_reduce_diagnol(lower, false);
	Matrix* z = matrix_unaugment(lower, lower->cols - 1);
	if (arr != NULL)
	{
		free(arr);
	}

	matrix_free(lower);

	matrix_augment(transpose, z);
	arr = matrix_reduce_lower(transpose, true, false);
	matrix_reduce_upper(transpose, false, false);
	matrix_reduce_diagnol(transpose, false);
	Matrix* sol = matrix_unaugment(transpose, transpose->cols - 1);
	if (arr != NULL)
	{
		free(arr);
	}

	matrix_free(transpose);

	printf("x:\n");
	matrix_print(sol, stdout);
	Matrix* mult = matrix_mult_matrix(m, sol);
	matrix_sub_matrix(b, mult);
	double error = matrix_get_magnitude(b);
	printf("Residual Norm: %lf\n", error);
	matrix_free(b);
	matrix_free(m);
	matrix_free(mult);
	matrix_free(sol);

	return 0;
}

int main(int argc, char* argv[])
{
	if (argc <= 1)
	{
		printf("usage: ./hw1 <problem number> <additional arguments>\n");
		return 1;
	}

	int problem_number = atoi(argv[1]);
	switch (problem_number)
	{
		default:
		{
			printf("invalid problem number\n");
			return 1;
			break;
		}
		case 2:
		{
			return problem_two(argc, argv);
			break;
		}
		case 3:
		{
			return problem_three(argc, argv);
			break;
		}
		case 4:
		{
			return problem_four(argc, argv);
			break;
		}
		case 5:
		{
			return problem_five(argc, argv);
			break;
		}
	}
}
