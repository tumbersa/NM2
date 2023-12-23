#pragma once
#include <fstream>

class Vector
{
	int size;
	double* values;
public:
	enum class CONDITIONALITY { good, bad };
	enum class DIMENSION { little, big };

	Vector();
	Vector(int _size);
	Vector(int _Size, int border1, int border2);
	Vector(int _Size, int border1, int border2, CONDITIONALITY con);
	Vector(int _Size, DIMENSION dim);
	Vector(std::ifstream& file);
	Vector(const Vector& a);
	~Vector();

	Vector& operator=(const Vector& a);
	Vector operator+(const Vector& a) const;
	Vector operator-(const Vector& a) const;
	Vector operator*(double scalar) const;
	double operator*(const Vector& a) const;

	void setValues(int pos, double val);
	double getValues(int pos) const;
	int getSize() const;
	double getNorm() const;
	void print() const;
};