/*!
* file: Matrix_Solution.cpp
* synopsis: Matrix_Methods test drive
* variant: 15
* created by Mykola Ozerov KV-53
*/
#include <fstream>
#include "MatrixMethods.h"

int main()
{
	double eps = 1e-8;
	MatrixMethods MyMethods;
	ifstream fin("input.txt");
	int N;

	if (!fin.is_open())
	{
		printf("Can't open file\n");
		return 0;
	}

	fin >> N;

	vector<vector<double>> A(N, vector<double>(N + 1));
	vector<double> X(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j <= N; j++)
			fin >> A[i][j];
	
	// Diagonally dominant matrix
	// Row 1: (1) - 3 * (4)
	// Row 4: 2 * (4) + (1) - (3)
	printf("Gaussian Elimination:\n\n");
	if (MyMethods.GaussianElimination(A, X))
		for (int i = 0; i < N; i++)
			printf("X%d = %.2lf\t", i + 1, X[i]);
	else
		printf("Can't find solution for Gaussian elimination");
	printf("\n");

	MyMethods.DirectIteration(A, X, eps);
	printf("\nDirect Iteration:\n\n");
	for (int i = 0; i < N; i++)
		printf("X%d = %.2lf\t", i + 1, X[i]);
	printf("\n");

	return 0;
}