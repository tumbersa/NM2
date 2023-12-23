#include <iostream>
#include <iomanip>
#include <Windows.h>
#include "Vector.h"
#include "Matrix2Gauss.h"

void answerCheck2();
void answerCheck3();

int main() {
	srand(time(0));
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	int SIZE = 8;
	Matrix A(SIZE, -10, 10);
	std::cout << "�������� �������:\n";
	A.print(); std::cout << std::endl;

	Vector X_exact(SIZE, 10, 20);
	std::cout << "������ X:\n";
	X_exact.print();
	Vector F(SIZE);
	F = A * X_exact;
	std::cout << "������ F:\n";
	F.print();

	Vector X(SIZE);
	X = A.solution(F, true,X_exact,true);
	std::cout << "����������� ������ X: ";
	X.print();
	std::cout << "����������� ������������ � ���������� �������: " << (X_exact - X).getNorm() << "\n";

	answerCheck2();
	answerCheck3();
}

void answerCheck2()
{
	std::cout << "\n" << std::setw(4) << "SIZE" << "   " << std::setw(3) << "k" << "   " << std::setw(4) << "�����������\n";
	int SIZE = 16;
	while (SIZE <= 1024)
	{
		for (int k = 2; k <= 6; k += 2)
		{
			Matrix A(SIZE, -10, 10, Vector::CONDITIONALITY::bad, k);
			Vector X_exact(SIZE, 10, 20);
			Vector F(SIZE);
			F = A * X_exact;

			Vector X(SIZE);
			X = A.solution(F, false);

			if (k == 2) std::cout << std::setw(4) << SIZE << "   ";
			else std::cout << std::setw(4) << " " << "   ";

			std::cout << std::setw(3) << k << "   ";
			std::cout << (X_exact - X).getNorm() << "\n";
		}

		SIZE *= 4;
	}
}
void answerCheck3()
{
	std::cout << "\n" << std::setw(4) << "SIZE" << "   " << "�����������\n";
	int SIZE = 16;
	int L = 3;
	while (SIZE <= 1024)
	{
		std::cout << std::setw(4) << SIZE << "   ";
		Matrix A(SIZE, 10, 100, Vector::CONDITIONALITY::good, 2);
		Vector X_exact(SIZE, 10, 20);
		Vector F(SIZE);
		F = A * X_exact;

		Vector X(SIZE);
		X = A.solution(F, false);
		std::cout << (X_exact - X).getNorm() << "\n";

		SIZE *= 4;
	}
}