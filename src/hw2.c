#include <stdio.h>
#include <string.h>
#include <matrix.h>

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

static int problem_two(int argc, char* argv[])
{
	return 0;
}

static int problem_three(int argc, char* argv[])
{
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
			return problem_two(argc, argv);
			break;
		}
		case 3:
		{
			return problem_three(argc, argv);
			break;
		}
	}
}
