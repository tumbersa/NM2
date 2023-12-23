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
			if (j >= i) // рандомятся элементы главной диагонали и над главной диагональю
				setDataValue(i, j, border1 + static_cast<double>(rand() % (border2 - border1 + 1)));
			else // тк матрица симметричная, берутся элементы над главной диагональю
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
			if (j > i) // рандомятся элементы над главной диагональю
				setDataValue(i, j, border1 + static_cast<double>(rand() % (border2 - border1 + 1)));
			else if (j == i) // находим элементы на главной диагонали, учитывая обусловленность
			{
				if (con == Vector::CONDITIONALITY::good)
					setDataValue(i, i, (border1 + static_cast<double>(rand() % (border2 - border1 + 1))) * pow(10, k));
				else
					setDataValue(i, i, (border1 + static_cast<double>(rand() % (border2 - border1 + 1))) / pow(10, k));
			}
			else //тк матрица симметричная, берутся элементы над главной диагональю
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
		std::cout << "Разные размерности. Умножение матрицы на вектор невозможно.\n";
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
	// вспомогательный массив, хранящий последовательность ведущих строк

	Vector P(size);

	// функция возвращает true - если номер строки(strNum) был уже задействован до текущего столбца(jNum), иначе false

	auto isStrUsed = [](const Vector& P, int strNum, int jNum) -> bool
	{
		for (int i = 0; i < jNum; i++)
			if (P.getValues(i) == strNum)
				return true;

		return false;
	};

	// функция возвращает номер строки с самым большим элементом в столбце jNum матрицы ourMatrix

	auto strWithMaxElem = [isStrUsed](const Vector& P, Matrix& ourMatrix, int jNum) -> int
	{
		int answer = 0;
		int maxEl = 0;

		for (int strNum = 0; strNum < ourMatrix.getSize(); strNum++)
		{
			if (!isStrUsed(P, strNum, jNum) && abs(ourMatrix.getDataValue(strNum, jNum)) >= maxEl) // если строка не была заюзана и элемент в этой строке больше maxEl, то
			{
				maxEl = abs(ourMatrix.getDataValue(strNum, jNum));
				answer = strNum;
			}
		}

		return answer;
	};


	// алгоритм решения системы методом Гаусса с частичной стратегией выбора ведущих элементов по столбцу


	for (int jNum = 0; jNum < size; jNum++)					// идём по каждой строке
	{
		int leadingLine = strWithMaxElem(P, *this, jNum);	// номер строки с ведущим элементом
		P.setValues(jNum, leadingLine);						// запоминаем этот номер на позиции jNum в P

		// приводим leadingLine к виду с единицей на главной диагонали

		double R = 1 / getDataValue(leadingLine, jNum);						// коэффицент деления
		setDataValue(leadingLine, jNum, 1);
		F.setValues(leadingLine, F.getValues(leadingLine) * R);				// умножаем зачение элемента F на коэфицент деления
		for (int j = jNum + 1; j < size; j++)
			setDataValue(leadingLine, j, getDataValue(leadingLine, j) * R);	// каждый элемент после jNum строки leadingLine умножается на коэффицент деления

		// обнуляем столбец jNum кроме элементов уже задействованных строк(хранятся в P)

		for (int strNum = 0; strNum < size; strNum++)
		{
			if (!isStrUsed(P, strNum, jNum + 1))
			{
				R = getDataValue(strNum, jNum);															// коэффицент который будем обнулять
				setDataValue(strNum, jNum, 0);
				F.setValues(strNum, F.getValues(strNum) - R * F.getValues(leadingLine));				// изменяем вектор F

				for (int j = jNum + 1; j < size; j++)
					setDataValue(strNum, j, getDataValue(strNum, j) - R * getDataValue(leadingLine, j));// изменяем каждый элемент после jNum строки strNum
			}
			
		}

		// необязательная печать

		if (print)
		{
			std::cout << "Матрица после " << jNum + 1 << "-ого шага:\n";
			adjointPrint(F);
		}

		if (printNorm) {
			std::cout << "Iteration " << jNum + 1 << ": \n";
			std::cout << (F - *this * X_beg).getNorm() << '\n';
		}
	}

	// необязательная печать

	if (print)
	{
		std::cout << "Строки с ведущими элементами: ";
		P.print();
	}

	// нахождение вектора решения

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