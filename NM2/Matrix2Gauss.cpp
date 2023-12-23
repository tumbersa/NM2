#include <iostream>
#include <iomanip>
#include "Matrix2Gauss.h"

Matrix::Matrix(int _size)
{
	size = _size;
	data = new double[size * size];
	for (int i = 0; i < size * size; i++) data[i] = 0;
}

Matrix::Matrix(int _size, int border1, int border2)
{
	size = _size;
	data = new double[size * size];

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (j >= i) // ���������� �������� ������� ��������� � ��� ������� ����������
				setDataValue(i, j, border1 + static_cast<double>(rand() % (border2 - border1 + 1)));
			else // �� ������� ������������, ������� �������� ��� ������� ����������
				setDataValue(i, j, getDataValue(j, i));
		}
	}
}

Matrix::Matrix(int _size, int border1, int border2, Vector::CONDITIONALITY con, int k)
{
	size = _size;
	data = new double[size * size];

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (j > i) // ���������� �������� ��� ������� ����������
				setDataValue(i, j, border1 + static_cast<double>(rand() % (border2 - border1 + 1)));
			else if (j == i) // ������� �������� �� ������� ���������, �������� ���������������
			{
				if (con == Vector::CONDITIONALITY::good)
					setDataValue(i, i, (border1 + static_cast<double>(rand() % (border2 - border1 + 1))) * pow(10, k));
				else
					setDataValue(i, i, (border1 + static_cast<double>(rand() % (border2 - border1 + 1))) / pow(10, k));
			}
			else //�� ������� ������������, ������� �������� ��� ������� ����������
				setDataValue(i, j, getDataValue(j, i));
		}
	}
}

Matrix::Matrix(const Matrix& x)
{
	size = x.size;
	data = new double[size * size];
	for (int i = 0; i < size * size; i++)
		data[i] = x.data[i];
}

Matrix::~Matrix()
{
	delete[] data;
}

Matrix& Matrix::operator=(const Matrix& x)
{
	if (&x != this)
	{
		size = x.size;
		data = new double[size * size];
		for (int i = 0; i < size * size; i++)
			data[i] = x.data[i];
	}
	return *this;
}

Vector Matrix::operator*(const Vector& x) const
{
	if (size != x.getSize())
	{
		std::cout << "������ �����������. ��������� ������� �� ������ ����������.\n";
		return x;
	}
	Vector answ(size);

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			answ.setValues(i, answ.getValues(i) + getDataValue(i, j) * x.getValues(j));

	return answ;
}

inline double Matrix::getDataValue(int i, int j) const
{
	return data[i * size + j];
}

inline void Matrix::setDataValue(int i, int j, double val)
{
	data[i * size + j] = val;
}

inline int Matrix::getSize() const
{
	return size;
}

void Matrix::print() const
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			std::cout << std::setw(8) << std::setprecision(4) << getDataValue(i, j);
		std::cout << "\n\n";
	}
}

void Matrix::adjointPrint(const Vector& F) const
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			std::cout << std::setw(8) << std::setprecision(4) << getDataValue(i, j);
		std::cout << std::setw(8) << "|" << F.getValues(i) << "\n\n";
	}
}

Vector Matrix::solution(Vector F, bool print, Vector X_beg, bool printNorm)
{
	// ��������������� ������, �������� ������������������ ������� �����

	Vector P(size);

	// ������� ���������� true - ���� ����� ������(strNum) ��� ��� ������������ �� �������� �������(jNum), ����� false

	auto isStrUsed = [](const Vector& P, int strNum, int jNum) -> bool
	{
		for (int i = 0; i < jNum; i++)
			if (P.getValues(i) == strNum)
				return true;

		return false;
	};

	// ������� ���������� ����� ������ � ����� ������� ��������� � ������� jNum ������� ourMatrix

	auto strWithMaxElem = [isStrUsed](const Vector& P, Matrix& ourMatrix, int jNum) -> int
	{
		int answer = 0;
		int maxEl = 0;

		for (int strNum = 0; strNum < ourMatrix.getSize(); strNum++)
		{
			if (!isStrUsed(P, strNum, jNum) && abs(ourMatrix.getDataValue(strNum, jNum)) >= maxEl) // ���� ������ �� ���� ������� � ������� � ���� ������ ������ maxEl, ��
			{
				maxEl = abs(ourMatrix.getDataValue(strNum, jNum));
				answer = strNum;
			}
		}

		return answer;
	};


	// �������� ������� ������� ������� ������ � ��������� ���������� ������ ������� ��������� �� �������


	for (int jNum = 0; jNum < size; jNum++)					// ��� �� ������ ������
	{
		int leadingLine = strWithMaxElem(P, *this, jNum);	// ����� ������ � ������� ���������
		P.setValues(jNum, leadingLine);						// ���������� ���� ����� �� ������� jNum � P

		// �������� leadingLine � ���� � �������� �� ������� ���������

		double R = 1 / getDataValue(leadingLine, jNum);						// ���������� �������
		setDataValue(leadingLine, jNum, 1);
		F.setValues(leadingLine, F.getValues(leadingLine) * R);				// �������� ������� �������� F �� ��������� �������
		for (int j = jNum + 1; j < size; j++)
			setDataValue(leadingLine, j, getDataValue(leadingLine, j) * R);	// ������ ������� ����� jNum ������ leadingLine ���������� �� ���������� �������

		// �������� ������� jNum ����� ��������� ��� ��������������� �����(�������� � P)

		for (int strNum = 0; strNum < size; strNum++)
		{
			if (!isStrUsed(P, strNum, jNum + 1))
			{
				R = getDataValue(strNum, jNum);															// ���������� ������� ����� ��������
				setDataValue(strNum, jNum, 0);
				F.setValues(strNum, F.getValues(strNum) - R * F.getValues(leadingLine));				// �������� ������ F

				for (int j = jNum + 1; j < size; j++)
					setDataValue(strNum, j, getDataValue(strNum, j) - R * getDataValue(leadingLine, j));// �������� ������ ������� ����� jNum ������ strNum
			}
			
		}

		// �������������� ������

		if (print)
		{
			std::cout << "������� ����� " << jNum + 1 << "-��� ����:\n";
			adjointPrint(F);
		}

		if (printNorm) {
			std::cout << "Iteration " << jNum + 1 << ": \n";
			std::cout << (F - *this * X_beg).getNorm() << '\n';
		}
	}

	// �������������� ������

	if (print)
	{
		std::cout << "������ � �������� ����������: ";
		P.print();
	}

	// ���������� ������� �������

	Vector X(size);
	for (int j = size - 1; j >= 0; j--)
	{
		double sum = 0;
		for (int jtemp = j + 1; jtemp < size; jtemp++)
			sum += getDataValue(P.getValues(j), jtemp) * X.getValues(jtemp);
		X.setValues(j, F.getValues(P.getValues(j)) - sum);
	}

	return X;
}