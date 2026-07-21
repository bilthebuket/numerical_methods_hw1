#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

#define TOLERANCE .00001

static int problem_one(int argc, char* argv[])
{
	if (argc != 6 && argc != 5)
	{
		printf("usage: ./hw2 1 <method> <max number of iterations> <system to solve> <omega (only required for SOR method)>\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[4]);
	if (m == NULL)
	{
		printf("Could not parse matrix\n");
		return 1;
	}

	int maxit = atoi(argv[3]);
	double tolerance = TOLERANCE;

	IterationHistory* ih;
	if (!strcmp(argv[2], "Jacobi"))
	{
		ih = stationary_solver(m, tolerance, maxit, JACOBI, 1.0);
	}
	else if (!strcmp(argv[2], "Gauss-Seidel"))
	{
		ih = stationary_solver(m, tolerance, maxit, GAUSS_SEIDEL, 1.0);
	}
	else if (!strcmp(argv[2], "SOR"))
	{
		double omega;
		sscanf(argv[5], "%lf", &omega);
		ih = stationary_solver(m, tolerance, maxit, SOR, omega);
	}
	else
	{
		printf("invalid solving method\n");
		matrix_free(m);
		return 1;
	}

	matrix_free(m);

	if (ih == NULL)
	{
		printf("solver failed\n");
		return 1;
	}

	for (int i = 0; i < 5; i++)
	{
		matrix_print(&(ih->xhist[i]), stdout);
		printf("\n");
	}
	matrix_print(&(ih->xhist[ih->iterations_executed - 1]), stdout);
	printf("\n");
	printf("Residual Norm: %lf\n", ih->rhist[ih->iterations_executed - 1]);
	printf("Iteration Count: %d\n", ih->iterations_executed);

	ih_free(ih);

	return 0;
}

static int problem_one_part_two(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("usage: ./hw2 2 <system to solve>\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}

	FILE* f = fopen("output.csv", "w");
	if (f == NULL)
	{
		printf("could not open file\n");
		matrix_free(m);
		return 1;
	}

	IterationHistory* ih1 = stationary_solver(m, TOLERANCE, 100, JACOBI, 1.0);
	IterationHistory* ih2 = stationary_solver(m, TOLERANCE, 100, GAUSS_SEIDEL, 1.0);
	IterationHistory* ih3 = stationary_solver(m, TOLERANCE, 100, SOR, .8);

	int max = ih1->iterations_executed;
	if (ih2->iterations_executed > max)
	{
		max = ih2->iterations_executed;
	}
	if (ih3->iterations_executed > max)
	{
		max = ih3->iterations_executed;
	}

	for (int i = 0; i < max; i++)
	{
		double one;
		double two;
		double three;

		if (i < ih1->iterations_executed)
		{
			one = ih1->rhist[i];
		}
		else
		{
			one = ih1->rhist[ih1->iterations_executed - 1];
		}
		if (i < ih2->iterations_executed)
		{
			two = ih2->rhist[i];
		}
		else
		{
			two = ih2->rhist[ih2->iterations_executed - 1];
		}
		if (i < ih3->iterations_executed)
		{
			three = ih3->rhist[i];
		}
		else
		{
			three = ih3->rhist[ih3->iterations_executed - 1];
		}

		fprintf(f, "%d,%lf,%lf,%lf\n", i + 1, one, two, three);
	}

	matrix_free(m);
	ih_free(ih1);
	ih_free(ih2);
	ih_free(ih3);
	fclose(f);
	return 0;
}

#define MAXIT 100

static int problem_three(int argc, char* argv[])
{
	if (argc != 4 && argc != 3)
	{
		printf("usage: ./hw2 3 <matrix to find dominant eigenvalue for> <inital vector (optional)>\n");
		return 1;
	}

	Matrix* m = matrix_create(argv[2]);
	if (m == NULL)
	{
		printf("could not parse matrix\n");
		return 1;
	}

	Matrix* x;
	if (argc == 4)
	{
		x = matrix_create(argv[3]);
		if (x == NULL)
		{
			printf("could not parse matrix\n");
			matrix_free(m);
			return 1;
		}
	}
	else
	{
		char buf[m->rows * 2];
		for (int i = 0; i < (m->rows - 1) * 2; i += 2)
		{
			buf[i] = '1';
			buf[i + 1] = ';';
		}
		buf[m->rows * 2 - 2] = '1';
		buf[m->rows * 2 - 1] = '\0';

		x = matrix_create(buf);
	}
	double eigen_value;

	double x1[MAXIT];
	double x2[MAXIT];
	double x3[MAXIT];
	double mus[MAXIT];
	double diffs[MAXIT];
	double norms[MAXIT];

	char buf[6];
	buf[0] = '1';
	buf[1] = ';';
	buf[2] = '1';
	buf[3] = ';';
	buf[4] = '1';
	buf[5] = '\0';
	Matrix* true_vector = matrix_create(buf);
	double true_value = 4.0;

	int i = 0;

	for (; i < MAXIT; i++)
	{
		Matrix* y = matrix_mult_matrix(m, x);
		double mu = matrix_infinity_norm(y);

		x1[i] = x->vals[0];
		x2[i] = x->vals[1];
		if (m->cols > 2)
		{
			x3[i] = x->vals[2];
		}
		mus[i] = mu;

		diffs[i] = fabs(mu - true_value);

		Matrix* dupe = matrix_duplicate(x);
		matrix_sub_matrix(dupe, true_vector);
		norms[i] = matrix_infinity_norm(dupe);
		matrix_free(dupe);

		matrix_divide_constant(y, mu);
		if (matrix_equals(y, x))
		{
			matrix_free(x);
			x = y;
			eigen_value = mu;
			break;
		}
		matrix_mult_constant(x, -1);
		if (matrix_equals(x, y))
		{
			matrix_free(x);
			x = y;
			eigen_value = mu * -1;
			break;
		}
		matrix_free(x);
		x = y;
	}

	if (i < MAXIT)
	{
		x1[i] = x->vals[0];
		x2[i] = x->vals[1];
		if (m->cols > 2)
		{
			x3[i] = x->vals[2];
		}
		mus[i] = eigen_value;

		diffs[i] = fabs(eigen_value - true_value);

		Matrix* dupe = matrix_duplicate(x);
		matrix_sub_matrix(dupe, true_vector);
		norms[i] = matrix_infinity_norm(dupe);
		matrix_free(dupe);
	}

	printf("eigen vector:\n");
	matrix_print(x, stdout);
	printf("eigen value: %lf\n", eigen_value);

	FILE* f = fopen("vectors.txt", "w");
	for (int j = 0; j < i; j++)
	{
		fprintf(f, "%lf,", x1[j]);
	}
	fprintf(f, "%lf;", x1[i]);

	for (int j = 0; j < i; j++)
	{
		fprintf(f, "%lf,", x2[j]);
	}
	fprintf(f, "%lf;", x2[i]);

	if (m->cols > 2)
	{
		for (int j = 0; j < i; j++)
		{
			fprintf(f, "%lf,", x3[j]);
		}
		fprintf(f, "%lf\n", x3[i]);
	}

	fclose(f);
	f = fopen("values.txt", "w");

	for (int j = 0; j < i; j++)
	{
		fprintf(f, "%lf,", mus[j]);
	}
	fprintf(f, "%lf\n", mus[i]);
	
	fclose(f);
	f = fopen("plot1.txt", "w");

	for (int j = 0; j < i; j++)
	{
		fprintf(f, "%lf,", diffs[j]);
	}
	fprintf(f, "%lf\n", diffs[i]);

	fclose(f);
	f = fopen("plot2.txt", "w");
	for (int j = 0; j < i; j++)
	{
		fprintf(f, "%lf,", norms[j]);
	}
	fprintf(f, "%lf\n", norms[i]);


	matrix_free(m);
	matrix_free(x);

	return 0;
}

int main(int argc, char* argv[])
{
	if (argc <= 1)
	{
		printf("usage: ./hw2 <problem number> <additional arguments>\n");
		return 1;
	}

	int problem_number = atoi(argv[1]);
	switch (problem_number)
	{
		case 1:
		{
			return problem_one(argc, argv);
			break;
		}
		case 2:
		{
			return problem_one_part_two(argc, argv);
			break;
		}
		case 3:
		{
			return problem_three(argc, argv);
			break;
		}
	}
}
