#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("usage: ./test <matrix 1> <matrix 2> <good, decent or bad>\n");
		return 1;
	}
	FILE* f = fopen(argv[1], "r");
	if (f == NULL)
	{
		return 1;
	}

	fseek(f, 0, SEEK_END);
	size_t size = ftell(f);
	fseek(f, 0, SEEK_SET);

	char* buf = malloc(sizeof(char) * (size + 1));
	if (buf == NULL)
	{
		fclose(f);
		return 1;
	}
	size_t bytes_read = fread(buf, sizeof(char), size, f);
	buf[bytes_read] = '\0';
	fclose(f);

	Matrix* m1 = matrix_create(buf);
	free(buf);

	f = fopen(argv[2], "r");
	if (f == NULL)
	{
		return 1;
	}

	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);

	buf = malloc(sizeof(char) * (size + 1));
	if (buf == NULL)
	{
		fclose(f);
		return 1;
	}
	bytes_read = fread(buf, sizeof(char), size, f);
	buf[bytes_read] = '\0';
	fclose(f);

	Matrix* m2 = matrix_create(buf);
	free(buf);

	Matrix* m3;
	if (!strcmp(argv[3], "good"))
	{
		m3 = matrix_mult_matrix(m1, m2);
	}
	else if (!strcmp(argv[3], "bad"))
	{
		m3 = matrix_mult_matrix_bad(m1, m2);
	}
	else if (!strcmp(argv[3], "multi"))
	{
		m3 = matrix_mult_matrix_multithread(m1, m2);
	}
	else
	{
		m3 = matrix_mult_matrix_decent(m1, m2);
	}
	matrix_print(m1, stdout);
	printf("\n");
	matrix_print(m2, stdout);
	printf("\n");
	matrix_print(m3, stdout);
	matrix_free(m1);
	matrix_free(m2);
	matrix_free(m3);
}
