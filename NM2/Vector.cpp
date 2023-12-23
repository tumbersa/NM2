#include "Vector.h"
#include <iostream>
#include <string>

Vector::Vector()
{
	size = 0;
	values = nullptr;
}

Vector::Vector(int _size)
{
	size = _size;
	values = new double[size];
	for (int i = 0; i < size; i++)
		values[i] = 0;
}

Vector::Vector(int _Size, int border1, int border2)
{
	size = _Size;
	values = new double[size];
	for (int i = 0; i < size; i++)
		values[i] = border1 + rand() % (border2 - border1);
}

Vector::Vector(int _Size, int border1, int border2, CONDITIONALITY con)
{
	size = _Size;
	values = new double[size];
	if (con == CONDITIONALITY::good)
	{
		for (int i = 0; i < size; i++)
			values[i] = border1 + rand() % (border2 - border1);
	}
	else
	{
		for (int i = 0; i < size; i++)
			values[i] = border1 / 2 + rand() % (border2 - border1) / 2;
	}
}

Vector::Vector(int _Size, DIMENSION dim)
{
	size = _Size;
	values = new double[size];
	if (dim == DIMENSION::big)
	{
		for (int i = 0; i < size; i++)
			values[i] = rand() % 1000 + 1;
	}
	else
	{
		for (int i = 0; i < size; i++)
			values[i] = rand() % 10 + 1;
	}
}

Vector::Vector(std::ifstream& file)
{
	std::string info;
	file >> size;
	values = new double[size];
	for (int i = 0; i < size; i++)
		file >> values[i];
}

Vector::Vector(const Vector& a)
{
	size = a.size;
	values = new double[size];
	for (int i = 0; i < size; i++)
		values[i] = a.values[i];
}

Vector::~Vector()
{
	delete values;
}

Vector& Vector::operator=(const Vector& a)
{
	if (&a != this)
	{
		size = a.size;
		values = new double[size];
		for (int i = 0; i < size; i++)
			values[i] = a.values[i];
	}
	return *this;
}

Vector Vector::operator+(const Vector& a) const
{
	if (size != a.size)
	{
		std::cerr << "Векторы имеют разную размерность! Сложение невозможно\n";
		std::exit(1);
	}
	else
	{
		Vector tmp(size);
		for (int i = 0; i < size; i++)
			tmp.values[i] = values[i] + a.values[i];
		return tmp;
	}
}

Vector Vector::operator-(const Vector& a) const
{
	if (size != a.size)
	{
		std::cerr << "Векторы имеют разную размерность! Вычитание невозможно\n";
		std::exit(1);
	}
	else
	{
		Vector tmp(size);
		for (int i = 0; i < size; i++)
			tmp.values[i] = values[i] - a.values[i];
		return tmp;
	}
}

Vector Vector::operator*(double scalar) const
{
	Vector tmp(size);
	for (int i = 0; i < size; i++)
		tmp.values[i] = values[i] * scalar;
	return tmp;
}

double Vector::operator*(const Vector& a) const
{
	if (size != a.size)
	{
		std::cerr << "Векторы имеют разную размерность! Умножение невозможно\n";
		std::exit(1);
	}
	else
	{
		double sum = 0;
		for (int i = 0; i < size; i++)
			sum += values[i] * a.values[i];
		return sum;
	}
}


void Vector::setValues(int pos, double val)
{
	values[pos] = val;
}

double Vector::getValues(int pos) const
{
	return values[pos];
}

int Vector::getSize() const
{
	return size;
}

double Vector::getNorm() const
{
	if (!values)
	{
		std::cerr << "Вектор не определён! Норма не найдена\n";
		std::exit(1);
	}
	else
	{
		double max = abs(values[0]);
		for (int i = 1; i < size; i++)
		{
			if (abs(values[i]) > max)
				max = abs(values[i]);
		}
		return max;
	}
}

void Vector::print() const
{
	for (int i = 0; i < size; i++)
	{
		if (i == 0) std::cout << values[i];
		else std::cout << " " << values[i];
	}
	std::cout << std::endl;
}