#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "matrix.h"

IterationHistory* ih_create(int maxit)
{
	IterationHistory* ih = malloc(sizeof(IterationHistory));
	if (ih == NULL)
	{
		return NULL;
	}

	char* arrays = malloc(sizeof(double) * maxit * 2 + sizeof(Matrix) * maxit);
	if (arrays == NULL)
	{
		free(ih);
		return NULL;
	}
	ih->xhist = (Matrix*) &arrays[0];
	ih->rhist = (double*) &arrays[sizeof(Matrix) * maxit];
	ih->diffhist = (double*) &arrays[sizeof(Matrix) * maxit + sizeof(double) * maxit];
	ih->iterations_executed = 0;
	return ih;
}

void ih_free(IterationHistory* ih)
{
	if (ih == NULL)
	{
		return;
	}

	for (int i = 0; i < ih->iterations_executed; i++)
	{
		free(ih->xhist[i].vals);
	}
	free(ih->xhist);
	free(ih);
}

IterationHistory* stationary_solver(Matrix* m, double tolerance, int maxit, int method, double omega)
{
	if (m == NULL || maxit <= 0 || method < 0 || method > MAX_VALID_METHOD)
	{
		return NULL;
	}

	IterationHistory* r = ih_create(maxit);
	if (r == NULL)
	{
		return NULL;
	}

	char buf[m->rows * 2];
	for (int i = 0; i < 2 * (m->rows - 1); i += 2)
	{
		buf[i] = '0';
		buf[i + 1] = ';';
	}
	buf[2 * m->rows - 2] = '0';
	buf[2 * m->rows - 1] = '\0';
	Matrix* xold = matrix_create(buf);
	Matrix* b = matrix_unaugment(m, m->cols - 1);

	matrix_copy(&(r->xhist[0]), xold);
	r->rhist[0] = matrix_infinity_norm(b);
	r->diffhist[0] = (double) NAN;
	r->iterations_executed++;

	for (int i = 0; i < maxit - 1; i++)
	{
		Matrix* xnew = matrix_create(buf);
		for (int j = 0; j < m->rows; j++)
		{
			double sum = 0.0;
			for (int k = j + 1; k < m->cols; k++)
			{
				sum += m->vals[j * m->cols + k] * xold->vals[k];
			}
			if (method == JACOBI)
			{
				for (int k = 0; k < j; k++)
				{
					sum += m->vals[j * m->cols + k] * xold->vals[k];
				}
			}
			else
			{
				for (int k = 0; k < j; k++)
				{
					sum += m->vals[j * m->cols + k] * xnew->vals[k];
				}
			}

			sum = (b->vals[j] - sum);
			sum /= m->vals[j * m->cols + j];
			sum *= omega;
			sum += (1 - omega) * xold->vals[j];
			xnew->vals[j] = sum;
		}

		Matrix* mult = matrix_mult_matrix(m, xnew);
		// this is backwards but we are only using it to compute the infinity norm which takes the absolute value anyways
		// also this means that b does not get mutated and can be reused
		matrix_sub_matrix(mult, b);
		double residual_norm = matrix_infinity_norm(mult);
		r->rhist[i + 1] = residual_norm;
		matrix_free(mult);

		Matrix* xnew_dupe = matrix_duplicate(xnew);
		matrix_sub_matrix(xnew_dupe, xold);
		double numerator = matrix_infinity_norm(xnew_dupe);
		matrix_free(xnew_dupe);
		double denominator = matrix_infinity_norm(xold);
		if (denominator < 1.0)
		{
			denominator = 1.0;
		}
		r->diffhist[i + 1] = numerator / denominator;

		matrix_copy(&(r->xhist[i + 1]), xnew);
		free(xold);
		xold = xnew;
		r->iterations_executed++;

		if (r->diffhist[i + 1] < tolerance)
		{
			break;
		}
	}

	free(xold);
	matrix_augment(m, b);
	return r;
}

#define NUM_SIG_DIGITS 3

void truncate_num(double* d)
{
	if (d == NULL || *d == 0.0)
	{
		return;
	}

	double num_digits_left_of_decimal = floor(log10(fabs(*d))) + 1;
	double multiplier = pow(10, num_digits_left_of_decimal - NUM_SIG_DIGITS);
	*d = trunc(*d / multiplier) * multiplier;
}

Matrix* matrix_create(char* s)
{
	if (s == NULL)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}

	int rows = 1;
	int cols = -1;
	int this_cols = 1;
	for (int i = 0; s[i] != '\0'; i++)
	{
		if (s[i] == ',')
		{
			this_cols++;
		}
		else if (s[i] == ';')
		{
			if (cols == -1)
			{
				cols = this_cols;
			}
			else if (cols != this_cols)
			{
				free(r);
				return NULL;
			}

			rows++;
			this_cols = 1;
		}
	}

	if (cols == -1)
	{
		cols = this_cols;
	}
	else if (cols != this_cols)
	{
		free(r);
		return NULL;
	}

	r->rows = rows;
	r->cols = cols;
	r->vals = malloc(sizeof(double) * rows * cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	int store = 0;
	int num_vals_read = 0;
	for (int i = 0; s[i] != '\0'; i++)
	{
		if (s[i] == ',')
		{
			s[i] = '\0';
			sscanf(&s[store], "%lf", &r->vals[num_vals_read]);
			s[i] = ',';
			store = i + 1;
			num_vals_read++;
		}
		else if (s[i] == ';')
		{
			s[i] = '\0';
			sscanf(&s[store], "%lf", &r->vals[num_vals_read]);
			s[i] = ';';
			store = i + 1;
			num_vals_read++;
		}
	}
	sscanf(&s[store], "%lf", &r->vals[num_vals_read]);

	return r;
}

Matrix* matrix_duplicate(Matrix* m)
{
	if (m == NULL)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}

	r->rows = m->rows;
	r->cols = m->cols;
	r->vals = malloc(sizeof(double) * r->rows * r->cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	for (int i = 0; i < r->rows * r->cols; i++)
	{
		r->vals[i] = m->vals[i];
	}
	return r;
}

void matrix_free(Matrix* m)
{
	if (m != NULL)
	{
		free(m->vals);
		free(m);
	}
}

void matrix_copy(Matrix* copy_to, Matrix* copy_from)
{
	if (copy_to == NULL || copy_from == NULL)
	{
		return;
	}

	copy_to->vals = copy_from->vals;
	copy_to->rows = copy_from->rows;
	copy_to->cols = copy_from->cols;
}

void matrix_add_matrix(Matrix* add_to, Matrix* add_from)
{
	if (add_to == NULL || add_from == NULL || add_to->rows != add_from->rows || add_to->cols != add_from->cols)
	{
		return;
	}

	for (int i = 0; i < add_to->rows * add_to->cols; i++)
	{
		add_to->vals[i] += add_from->vals[i];
	}
}

void matrix_sub_matrix(Matrix* sub_from, Matrix* to_sub)
{
	if (sub_from == NULL || to_sub == NULL || sub_from->rows != to_sub->rows || sub_from->cols != to_sub->cols)
	{
		return;
	}

	for (int i = 0; i < sub_from->rows * to_sub->cols; i++)
	{
		sub_from->vals[i] -= to_sub->vals[i];
	}
}

#define CHUNK_SIZE 16

static void* matrix_mult_matrix_thread_func(void* v)
{
	void** vals = (void**) v;
	Matrix* left = (Matrix*) vals[0];
	Matrix* right = (Matrix*) vals[1];
	int* chunk_col_start = (int*) vals[2];
	int* chunk_col_end = (int*) vals[3];
	Matrix* r = (Matrix*) vals[4];

	for (int chunk_col = *chunk_col_start; chunk_col <= *chunk_col_end; chunk_col++)
	{
		for (int chunk_row = 0; chunk_row <= right->rows / CHUNK_SIZE; chunk_row++)
		{
			int num_cols = CHUNK_SIZE;
			int num_rows = CHUNK_SIZE;
			if (chunk_col == right->cols / CHUNK_SIZE)
			{
				num_cols = right->cols % CHUNK_SIZE;
			}
			if (chunk_row == right->rows / CHUNK_SIZE)
			{
				num_rows = right->rows % CHUNK_SIZE;
			}
			double active_chunk[num_cols * num_rows];
			for (int row_offset = 0; row_offset < num_rows; row_offset++)
			{
				for (int col_offset = 0; col_offset < num_cols; col_offset++)
				{
					active_chunk[col_offset * num_rows + row_offset] = right->vals[(chunk_row * CHUNK_SIZE + row_offset) * right->cols + chunk_col * CHUNK_SIZE + col_offset];
				}
			}

			for (int left_chunk_row = 0; left_chunk_row <= left->rows / CHUNK_SIZE; left_chunk_row++)
			{
				int left_num_rows = CHUNK_SIZE;
				if (left_chunk_row == left->rows / CHUNK_SIZE)
				{
					left_num_rows = left->rows % CHUNK_SIZE;
				}
				for (int left_row = 0; left_row < left_num_rows; left_row++)
				{
					for (int active_col = 0; active_col < num_cols; active_col++)
					{
						for (int offset = 0; offset < num_rows; offset++)
						{
							r->vals[(left_chunk_row * CHUNK_SIZE + left_row) * r->cols + chunk_col * CHUNK_SIZE + active_col] 
								+= left->vals[(left_chunk_row * CHUNK_SIZE + left_row) * left->cols + chunk_row * CHUNK_SIZE + offset] * active_chunk[active_col * num_rows + offset];
						}
					}
				}
			}
		}
	}

	return r;
}

#define NUM_THREADS 8

Matrix* matrix_mult_matrix_multithread(Matrix* left, Matrix* right)
{
	if (left == NULL || right == NULL || left->cols != right->rows)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->rows = left->rows;
	r->cols = right->cols;
	r->vals = calloc(r->rows * r->cols, sizeof(double));
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	pthread_t threads[NUM_THREADS];
	void* args[NUM_THREADS][5];
	int starts[NUM_THREADS];
	int ends[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++)
	{
		args[i][0] = left;
		args[i][1] = right;
		starts[i] = right->cols / CHUNK_SIZE / NUM_THREADS * i;
		ends[i] = right->cols / CHUNK_SIZE / NUM_THREADS * (i + 1) - 1;
		if (i == NUM_THREADS - 1)
		{
			ends[i] += (right->cols / CHUNK_SIZE) % NUM_THREADS;
			if (right->cols % CHUNK_SIZE != 0)
			{
				ends[i]++;
			}
		}
		args[i][2] = &starts[i];
		args[i][3] = &ends[i];
		args[i][4] = r;

		if (ends[i] < starts[i])
		{
			threads[i] = 0;
		}
		else
		{
			int result = pthread_create(&threads[i], NULL, &matrix_mult_matrix_thread_func, args[i]);
			if (result != 0)
			{
				for (int j = 0; j < i; j++)
				{
					pthread_join(threads[j], NULL);
					matrix_free(r);
				}
				return NULL;
			}
		}
	}
	for (int i = 0; i < NUM_THREADS; i++)
	{
		if (threads[i] == 0)
		{
			continue;
		}
		pthread_join(threads[i], NULL);
	}
	return r;
}

Matrix* matrix_mult_matrix(Matrix* left, Matrix* right)
{
	if (left == NULL || right == NULL || left->cols != right->rows)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}

	r->rows = left->rows;
	r->cols = right->cols;
	r->vals = calloc(r->rows * r->cols, sizeof(double));
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	for (int left_chunk_row = 0; left_chunk_row <= left->rows / CHUNK_SIZE; left_chunk_row++)
	{
		for (int chunk_col = 0; chunk_col <= right->cols / CHUNK_SIZE; chunk_col++)
		{
			for (int chunk_row = 0; chunk_row <= right->rows / CHUNK_SIZE; chunk_row++)
			{
				int num_cols = CHUNK_SIZE;
				int num_rows = CHUNK_SIZE;
				if (chunk_col == right->cols / CHUNK_SIZE)
				{
					num_cols = right->cols % CHUNK_SIZE;
				}
				if (chunk_row == right->rows / CHUNK_SIZE)
				{
					num_rows = right->rows % CHUNK_SIZE;
				}
				double active_chunk[num_cols * num_rows];
				for (int row_offset = 0; row_offset < num_rows; row_offset++)
				{
					for (int col_offset = 0; col_offset < num_cols; col_offset++)
					{
						active_chunk[col_offset * num_rows + row_offset] = right->vals[(chunk_row * CHUNK_SIZE + row_offset) * right->cols + chunk_col * CHUNK_SIZE + col_offset];
					}
				}
				int left_num_rows = CHUNK_SIZE;
				if (left_chunk_row == left->rows / CHUNK_SIZE)
				{
					left_num_rows = left->rows % CHUNK_SIZE;
				}
				for (int left_row = 0; left_row < left_num_rows; left_row++)
				{
					for (int active_col = 0; active_col < num_cols; active_col++)
					{
						for (int offset = 0; offset < num_rows; offset++)
						{
							r->vals[(left_chunk_row * CHUNK_SIZE + left_row) * r->cols + chunk_col * CHUNK_SIZE + active_col] 
								+= left->vals[(left_chunk_row * CHUNK_SIZE + left_row) * left->cols + chunk_row * CHUNK_SIZE + offset] * active_chunk[active_col * num_rows + offset];
						}
					}
				}
			}
		}
	}

	return r;
}

Matrix* matrix_mult_matrix_decent(Matrix* left, Matrix* right)
{
	if (left == NULL || right == NULL || left->cols != right->rows)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}

	r->rows = left->rows;
	r->cols = right->cols;
	r->vals = calloc(r->rows * r->cols, sizeof(double));
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	double active_column[CHUNK_SIZE];

	int i = 0;
	for (; i < right->rows / CHUNK_SIZE; i++)
	{
		int j = 0;
		for (; j < right->cols / CHUNK_SIZE; j++)
		{
			for (int k = 0; k < CHUNK_SIZE; k++)
			{
				for (int w = 0; w < CHUNK_SIZE; w++)
				{
					active_column[w] = right->vals[(i * CHUNK_SIZE + w) * right->cols + j * CHUNK_SIZE + k];
				}

				for (int w = 0; w < left->rows; w++)
				{
					for (int h = 0; h < CHUNK_SIZE; h++)
					{
						r->vals[w * r->cols + j * CHUNK_SIZE + k] += left->vals[w * left->cols + i * CHUNK_SIZE + h] * active_column[h];
					}
				}
			}
		}
		for (int k = 0; k < right->cols % CHUNK_SIZE; k++)
		{
			for (int w = 0; w < CHUNK_SIZE; w++)
			{
				active_column[w] = right->vals[(i * CHUNK_SIZE + w) * right->cols + j * CHUNK_SIZE + k];
			}

			for (int w = 0; w < left->rows; w++)
			{
				for (int h = 0; h < CHUNK_SIZE; h++)
				{
					r->vals[w * r->cols + j * CHUNK_SIZE + k] += left->vals[w * left->cols + i * CHUNK_SIZE + h] * active_column[h];
				}
			}
		}
	}
	int j = 0;
	for (; j < right->cols / CHUNK_SIZE; j++)
	{
		for (int k = 0; k < CHUNK_SIZE; k++)
		{
			for (int w = 0; w < right->rows % CHUNK_SIZE; w++)
			{
				active_column[w] = right->vals[(i * CHUNK_SIZE + w) * right->cols + j * CHUNK_SIZE + k];
			}

			for (int w = 0; w < left->rows; w++)
			{
				for (int h = 0; h < right->rows % CHUNK_SIZE; h++)
				{
					r->vals[w * r->cols + j * CHUNK_SIZE + k] += left->vals[w * left->cols + i * CHUNK_SIZE + h] * active_column[h];
				}
			}
		}
	}
	for (int k = 0; k < right->cols % CHUNK_SIZE; k++)
	{
		for (int w = 0; w < right->rows % CHUNK_SIZE; w++)
		{
			active_column[w] = right->vals[(i * CHUNK_SIZE + w) * right->cols + j * CHUNK_SIZE + k];
		}

		for (int w = 0; w < left->rows; w++)
		{
			for (int h = 0; h < right->rows % CHUNK_SIZE; h++)
			{
				r->vals[w * r->cols + j * CHUNK_SIZE + k] += left->vals[w * left->cols + i * CHUNK_SIZE + h] * active_column[h];
			}
		}
	}

	return r;
}

Matrix* matrix_mult_matrix_bad(Matrix* left, Matrix* right)
{
	if (left == NULL || right == NULL || left->cols != right->rows)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}

	r->rows = left->rows;
	r->cols = right->cols;
	r->vals = calloc(r->rows * r->cols, sizeof(double));
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	for (int i = 0; i < right->cols; i++)
	{
		for (int j = 0; j < left->rows; j++)
		{
			for (int k = 0; k < left->cols; k++)
			{
				r->vals[j * r->cols + i] += left->vals[j * left->cols + k] * right->vals[k * right->cols + i];
			}
		}
	}

	return r;
}

void matrix_print(Matrix* m, FILE* f)
{
	if (m == NULL || f == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j < m->cols; j++)
		{
			fprintf(f, "%lf ", m->vals[i * m->cols + j]);
		}
		fprintf(f, "\n");
	}
}

void matrix_add_constant(Matrix* m, double d)
{
	if (m == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows * m->cols; i++)
	{
		m->vals[i] += d;
	}
}

void matrix_sub_constant(Matrix* m, double d)
{
	if (m == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows * m->cols; i++)
	{
		m->vals[i] -= d;
	}
}

void matrix_mult_constant(Matrix* m, double d)
{
	if (m == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows * m->cols; i++)
	{
		m->vals[i] *= d;
	}
}

void matrix_divide_constant(Matrix* m, double d)
{
	if (m == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows * m->cols; i++)
	{
		m->vals[i] /= d;
	}
}

void matrix_swap_row(Matrix* m, int row1, int row2)
{
	if (m == NULL || row1 < 0 || row2 < 0 || row1 >= m->rows || row2 >= m->rows || row1 == row2)
	{
		return;
	}

	for (int i = 0; i < m->cols; i++)
	{
		double store = m->vals[row1 * m->cols + i];
		m->vals[row1 * m->cols + i] = m->vals[row2 * m->cols + i];
		m->vals[row2 * m->cols + i] = store;
	}
}

void matrix_swap_rows(Matrix* m, int* swaps)
{
	if (m == NULL || swaps == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows; i++)
	{
		if (swaps[i] != i)
		{
			for (int j = 0; j < m->cols; j++)
			{
				double store = m->vals[swaps[i] * m->cols + j];
				m->vals[swaps[i] * m->cols + j] = m->vals[swaps[swaps[i]] * m->cols + j];
				m->vals[swaps[swaps[i]] * m->cols + j] = store;
			}
			int store = swaps[i];
			swaps[i] = swaps[swaps[i]];
			swaps[store] = store;
		}
	}
}

void matrix_add_row(Matrix* m, int to_apply, int apply_to, double multiplier, bool truncate)
{
	if (m == NULL || to_apply < 0 || apply_to < 0 || to_apply >= m->rows || apply_to >= m->rows)
	{
		return;
	}

	for (int i = 0; i < m->cols; i++)
	{
		double multiplied = m->vals[to_apply * m->cols + i] * multiplier;
		if (truncate)
		{
			truncate_num(&multiplied);
		}
		m->vals[apply_to * m->cols + i] += multiplied;
		if (truncate)
		{
			truncate_num(&m->vals[apply_to * m->cols + i]);
		}
	}
}

void matrix_mult_row(Matrix* m, int row, double val, bool truncate)
{
	if (m == NULL || row < 0 || row >= m->rows)
	{
		return;
	}

	for (int i = 0; i < m->cols; i++)
	{
		m->vals[row * m->cols + i] *= val;
		if (truncate)
		{
			truncate_num(&m->vals[row * m->cols + i]);
		}
	}
}

void matrix_divide_row(Matrix*m, int row, double val, bool truncate)
{
	if (m == NULL || row < 0 || row >= m->rows)
	{
		return;
	}

	for (int i = 0; i < m->cols; i++)
	{
		m->vals[row * m->cols + i] /= val;
		if (truncate)
		{
			truncate_num(&m->vals[row * m->cols + i]);
		}
	}
}

void matrix_reduce_diagnol(Matrix* m, bool truncate)
{
	if (m == NULL)
	{
		return;
	}

	int min = m->cols;
	if (m->rows < min)
	{
		min = m->rows;
	}

	for (int i = 0; i < min; i++)
	{
		matrix_divide_row(m, i, m->vals[i * m->cols + i], truncate);
	}
}

void matrix_reduce_upper(Matrix* m, bool truncate, bool reduce_diagnol)
{
	if (m == NULL)
	{
		return;
	}

	int min = m->cols;
	if (m->rows < min)
	{
		min = m->rows;
	}
	for (int i = min - 1; i >= 0; i--)
	{
		if (reduce_diagnol)
		{
			matrix_divide_row(m, i, m->vals[i * m->cols + i], truncate);
		}
		for (int j = i - 1; j >= 0; j--)
		{
			double factor = (m->vals[j * m->cols + i] / m->vals[i * m->cols + i]) * -1;
			if (truncate)
			{
				truncate_num(&factor);
			}
			matrix_add_row(m, i, j, factor, truncate);
		}
	}
}

int* matrix_reduce_lower(Matrix* m, bool partial_pivot, bool truncate)
{
	if (m == NULL)
	{
		return NULL;
	}

	int* r = malloc(sizeof(int) * m->rows);
	if (r == NULL)
	{
		return NULL;
	}
	for (int i = 0; i < m->rows; i++)
	{
		r[i] = i;
	}

	int min = m->cols;
	if (m->rows < m->cols)
	{
		min = m->rows;
	}
	for (int i = 0; i < min; i++)
	{
		if (partial_pivot || m->vals[i * m->cols + i] == 0.0)
		{
			double max = fabs(m->vals[i * m->cols + i]);
			int index = i;
			for (int j = i + 1; j < m->rows; j++)
			{
				if (fabs(m->vals[j * m->cols + i]) > max)
				{
					max = fabs(m->vals[j * m->cols + i]);
					index = j;
				}
			}
			if (max == 0.0)
			{
				continue;
			}
			matrix_swap_row(m, i, index);
			int store = r[i];
			r[i] = r[index];
			r[index] = store;
		}

		for (int j = i + 1; j < m->rows; j++)
		{
			double factor = (m->vals[j * m->cols + i] / m->vals[i * m->cols + i]) * -1;
			if (truncate)
			{
				truncate_num(&factor);
			}
			matrix_add_row(m, i, j, factor, truncate);
		}
	}

	return r;
}

int* matrix_decompose(Matrix* m, bool partial_pivot, bool truncate)
{
	if (m == NULL)
	{
		return NULL;
	}

	int* r = malloc(sizeof(int) * m->rows);
	if (r == NULL)
	{
		return NULL;
	}
	for (int i = 0; i < m->rows; i++)
	{
		r[i] = i;
	}

	int min = m->cols;
	if (m->rows < m->cols)
	{
		min = m->rows;
	}
	for (int i = 0; i < min; i++)
	{
		if (partial_pivot || m->vals[i * m->cols + i] == 0.0)
		{
			double max = fabs(m->vals[i * m->cols + i]);
			int index = i;
			for (int j = i + 1; j < m->rows; j++)
			{
				if (fabs(m->vals[j * m->cols + i]) > max)
				{
					max = fabs(m->vals[j * m->cols + i]);
					index = j;
				}
			}
			if (max == 0.0)
			{
				continue;
			}
			matrix_swap_row(m, i, index);
			int store = r[i];
			r[i] = r[index];
			r[index] = store;
		}

		for (int j = i + 1; j < m->rows; j++)
		{
			double factor = (m->vals[j * m->cols + i] / m->vals[i * m->cols + i]) * -1;
			if (truncate)
			{
				truncate_num(&factor);
			}
			m->vals[j * m->cols + i] /= m->vals[i * m->cols + i];
			for (int k = i + 1; k < m->cols; k++)
			{
				m->vals[j * m->cols + k] += m->vals[i * m->cols + k] * factor;
				if (truncate)
				{
					truncate_num(&m->vals[j * m->cols + k]);
				}
			}
		}
	}

	return r;
}

void matrix_cholesky_decompose(Matrix* m)
{
	if (m == NULL || m->rows != m->cols)
	{
		return;
	}

	for (int i = 0; i < m->rows; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += pow(m->vals[i * m->cols + j], 2);
		}
		m->vals[i * m->cols + i] = pow(m->vals[i * m->cols + i] - sum, .5);
		for (int j = i + 1; j < m->rows; j++)
		{
			sum = 0.0;
			for (int k = 0; k < i; k++)
			{
				sum += m->vals[i * m->cols + k] * m->vals[j * m->cols + k];
			}
			m->vals[j * m->cols + i] = (m->vals[i * m->cols + j] - sum) / m->vals[i * m->cols + i];
		}
	}
}

void matrix_invert(Matrix* m)
{
	if (m == NULL || m->rows != m->cols)
	{
		return;
	}

	Matrix* to_augment = matrix_make_identity(m->rows);
	if (to_augment == NULL)
	{
		return;
	}

	int store = m->cols;
	matrix_augment(m, to_augment);
	matrix_free(to_augment);
	int* arr = matrix_reduce_lower(m, true, false);
	if (arr != NULL)
	{
		matrix_swap_rows(m, arr);
		free(arr);
	}
	Matrix* inverted = matrix_unaugment(m, store);
	if (inverted == NULL)
	{
		return;
	}
	double* store2 = m->vals;
	m->vals = inverted->vals;
	inverted->vals = store2;
	matrix_free(inverted);
}

void matrix_transpose(Matrix* m)
{
	if (m == NULL)
	{
		return;
	}
	
	double* new_vals = malloc(sizeof(double) * m->rows * m->cols);
	if (new_vals == NULL)
	{
		return;
	}

	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j < m->cols; j++)
		{
			new_vals[j * m->rows + i] = m->vals[i * m->cols + j];
		}
	}
	free(m->vals);
	m->vals = new_vals;
	int store = m->cols;
	m->cols = m->rows;
	m->rows = store;
}

Matrix* matrix_get_lower(Matrix* m)
{
	if (m == NULL)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->vals = malloc(sizeof(double) * m->rows * m->cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}
	r->cols = m->cols;
	r->rows = m->rows;

	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j < i; j++)
		{
			r->vals[i * m->cols + j] = m->vals[i * m->cols + j];
		}
		for (int j = i; j < m->cols; j++)
		{
			r->vals[i * m->cols + j] = 0.0;
		}
	}

	return r;
}

Matrix* matrix_get_upper(Matrix* m)
{
	if (m == NULL)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->vals = malloc(sizeof(double) * m->rows * m->cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}
	r->cols = m->cols;
	r->rows = m->rows;

	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			r->vals[i * m->cols + j] = 0.0;
		}
		for (int j = i + 1; j < m->cols; j++)
		{
			r->vals[i * m->cols + j] = m->vals[i * m->cols + j];
		}
	}

	return r;
}

Matrix* matrix_get_diagnol(Matrix* m)
{
	if (m == NULL)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->vals = malloc(sizeof(double) * m->rows * m->cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}
	r->cols = m->cols;
	r->rows = m->rows;

	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j < i; j++)
		{
			r->vals[i * m->cols + j] = 0.0;
		}
		r->vals[i * m->cols + i] = m->vals[i * m->cols + i];
		for (int j = i + 1; j < m->cols; j++)
		{
			r->vals[i * m->cols + j] = 0.0;
		}
	}

	return r;
}

Matrix* matrix_make_identity(int size)
{
	if (size <= 0)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->vals = calloc(size * size, sizeof(double));
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}
	r->rows = size;
	r->cols = size;

	for (int i = 0; i < size; i++)
	{
		r->vals[i * size + i] = 1.0;
	}

	return r;
}

void matrix_augment(Matrix* add_to, Matrix* add_from)
{
	if (add_to == NULL || add_from == NULL || add_to->rows != add_from->rows)
	{
		return;
	}

	double* new_vals = malloc(sizeof(double) * (add_to->rows * add_to->cols + add_from->rows * add_from->cols));
	if (new_vals == NULL)
	{
		return;
	}

	for (int i = 0; i < add_to->rows; i++)
	{
		for (int j = 0; j < add_to->cols; j++)
		{
			new_vals[i * (add_to->cols + add_from->cols) + j] = add_to->vals[i * add_to->cols + j];
		}
		for (int j = 0; j < add_from->cols; j++)
		{
			new_vals[i * (add_to->cols + add_from->cols) + j + add_to->cols] = add_from->vals[i * add_from->cols + j];
		}
	}

	free(add_to->vals);
	add_to->vals = new_vals;
	add_to->cols += add_from->cols;
	matrix_free(add_from);
}

Matrix* matrix_unaugment(Matrix* m, int cutoff_column)
{
	if (m == NULL || cutoff_column <= 0 || cutoff_column >= m->cols)
	{
		return NULL;
	}

	Matrix* r = malloc(sizeof(Matrix));
	if (r == NULL)
	{
		return NULL;
	}
	r->rows = m->rows;
	r->cols = m->cols - cutoff_column;
	r->vals = malloc(sizeof(double) * r->rows * r->cols);
	if (r->vals == NULL)
	{
		free(r);
		return NULL;
	}

	for (int i = 0; i < r->rows; i++)
	{
		for (int j = 0; j < r->cols; j++)
		{
			r->vals[i * r->cols + j] = m->vals[i * m->cols + j + cutoff_column];
		}
	}

	double* new_vals = malloc(sizeof(double) * m->rows * cutoff_column);
	if (new_vals == NULL)
	{
		return r;
	}
	for (int i = 0; i < m->rows; i++)
	{
		for (int j = 0; j < cutoff_column; j++)
		{
			new_vals[i * cutoff_column + j] = m->vals[i * m->cols + j];
		}
	}
	free(m->vals);
	m->vals = new_vals;
	m->cols = cutoff_column;

	return r;
}

double matrix_get_magnitude(Matrix* m)
{
	if (m == NULL)
	{
		return (double) NAN;
	}

	double r = 0.0;

	for (int i = 0; i < m->rows * m->cols; i++)
	{
		r += pow(m->vals[i], 2);
	}

	return pow(r, .5);
}

double matrix_infinity_norm(Matrix* m)
{
	if (m == NULL)
	{
		return (double) NAN;
	}

	double r = 0.0;

	for (int i = 0; i < m->rows; i++)
	{
		double this_row = 0.0;

		for (int j = 0; j < m->cols; j++)
		{
			this_row += fabs(m->vals[i * m->cols + j]);
		}

		if (this_row > r)
		{
			r = this_row;
		}
	}

	return r;
}

#define EQUALITY_THRESHOLD .000001

bool matrix_equals(Matrix* m1, Matrix* m2)
{
	if (m1 == NULL || m2 == NULL)
	{
		return false;
	}

	if (m1 == m2)
	{
		return true;
	}

	if (m1->cols != m2->cols || m1->rows != m2->rows)
	{
		return false;
	}

	for (int i = 0; i < m1->rows; i++)
	{
		for (int j = 0; j < m1->cols; j++)
		{
			if (fabs(m1->vals[i * m1->cols + j] - m2->vals[i * m1->cols + j]) > EQUALITY_THRESHOLD)
			{
				return false;
			}
		}
	}

	return true;
}
