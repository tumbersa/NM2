#pragma once
#include "Vector.h"

class Matrix
{
	int size;
	double* data;

public:
	Matrix(int _size);
	Matrix(int _size, int border1, int border2);
	Matrix(int _size, int border1, int border2, Vector::CONDITIONALITY con, int k);
	Matrix(const Matrix&);
	~Matrix();

	Matrix& operator=(const Matrix&);
	Vector operator*(const Vector&) const;

	inline double getDataValue(int i, int j) const;
	inline void setDataValue(int i, int j, double val);
	inline int getSize() const;

	void print() const;
	void adjointPrint(const Vector&) const;
	Vector solution(Vector F, bool print,Vector X_beg = Vector(),bool printNorm = false);
};