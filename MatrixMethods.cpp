/*!
* file: Matrix_Methods.cpp
* synopsis: Func definitions for Matrix_Methods
* variant: 15
* created by Mykola Ozerov KV-53
*/
#include "MatrixMethods.h"

bool MatrixMethods::GaussianElimination(vector<vector<double>> A, vector<double> &X)
{
	double tmp;
	int N = (int)X.size();

	//Direct elimination
	for (int i = 0; i < N; i++)
	{
		if (A[i][i] == 0)
			return false;

		tmp = A[i][i];
		for (int j = 0; j <= N; j++)
			A[i][j] /= tmp;
		
		for (int j = i + 1; j < N; j++)
		{
			tmp = A[j][i];
			for (int k = 0; k <= N; k++)
				A[j][k] = A[j][k] - tmp * A[i][k];
		}
	}

	//Backward substitution
	for (int i = N - 1; i >= 0; i--)
	{
		X[i] = A[i][N];
		for (int j = i + 1; j < N; j++)
			X[i] -= A[i][j] * X[j];
	}
	return true;
}

void MatrixMethods::DirectIteration(vector<vector<double>> A, vector<double> &X, double eps)
{
	double q, sum, tmp, fractionQ, norm;
	int N = (int)X.size();
	vector<double> prewX(N);

	//Conversion of matrix
	for (int i = 0; i < N; i++)
	{
		tmp = A[i][i];
		for (int j = 0; j < N; j++)
			if (j == i)
				A[i][j] /= tmp;
			else 
				A[i][j] /= -tmp;
		A[i][N] /= tmp;
	}
	
	//Initial approximation
	for (int i = 0; i < N; i++)
		prewX[i] = A[i][N];
	
	q = 0;
	for (int i = 0; i < N; i++)
	{		
		sum = 0;
		for (int j = 0; j < N; j++)
			if (i != j)
				sum += abs(A[i][j]);
		if (sum > q)
			q = sum;
	}
	fractionQ = (1 - q) / q;
	do
	{
		//Next		approximation
		for (int i = 0; i < N; i++)
		{	
			X[i] = A[i][N];
			for (int j = 0; j < N; j++)
				if (i != j)
					X[i] += A[i][j] * prewX[j];
		}

		//Norm calculation
		norm = 0;
		for (int i = 0; i < N; i++)
		{
			tmp = abs(X[i] - prewX[i]);
			if (tmp > norm)
				norm = tmp;
		}
		for (int i = 0; i < N; i++)
			prewX[i] = X[i];

	} while (norm > fractionQ * eps);
}