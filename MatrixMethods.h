/*!
* file: MatrixMethods.h
* synopsis: Common declarations for Matrix_Methods
* variant: 15
* created by Mykola Ozerov KV-53
*/
#pragma once
#include <vector>

using namespace std;

class MatrixMethods {
public:
	MatrixMethods() {};
	bool GaussianElimination(vector<vector<double>> A, vector<double> &X);
	void DirectIteration(vector<vector<double>> A, vector<double> &X, double eps);
	~MatrixMethods() {};
};