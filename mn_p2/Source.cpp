#include <iostream>
#include <math.h>
#include <ctime>
#include "../LinearEquationSolver/LinearEquationSolver.h"

#define N 977
#define A1 8
#define A2 -1
#define A3 -1
#define F 0
#define BLAD 1E-9

//#define DEBUG

#ifdef DEBUG
#define N 4
#endif

using namespace std;

void wypelnijMacierz(double* macierz)
{
#ifdef DEBUG
	macierz = {
		1.2, 2.6, -0.1, 1.5,
		4.5, 9.8, -0.4, 5.7,
		0.1, -0.1, -0.3, -3.5,
		4.5, -5.2, 4.2, -3.4,
	}

#else
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			macierz[i*N + j] = 0;
			if (i == j)
			{
				macierz[i*N + j] = A1;
				continue;
			}
			if (i + 1 == j || i == j + 1)
			{
				macierz[i*N + j] = A2;
				continue;
			}
			if (i + 2 == j || i == j + 2)
			{
				macierz[i*N + j] = A3;
			}
		}
	}
#endif
}

void wypelnijWektor(double* wektor)
{
#ifdef DEBUG
	wektor[0] = 13.15;
	wektor[1] = 49.84;
	wektor[2] = -14.08;
	wektor[3] = -46.51;
#else
	for (int i = 0; i < N; i++)
	{
		wektor[i] = sin(static_cast<float>(i + 1)*(F + 1) / 50);
	}
#endif
}

int main()
{
	auto macierz = new double[N*N];
	auto x = new double[N];
	auto wektor = new double[N];

	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);

	CLinearEquationSolver solver;

	clock_t begin = clock();
	solver.algorytmGaussa(macierz, wektor, x, N);
	cout << "Norma z residuum: " << solver.stopienBledu(x, macierz, wektor, N) << endl;
	clock_t end = clock();
	cout << "Czas: " << double(end - begin) / CLOCKS_PER_SEC * 1000 << endl;
	//cout << (double)(end - begin) << endl;
	//print(x);

	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);
	begin = clock();
	int iteracje = solver.algorytmJacobiego(macierz, wektor, x, N, BLAD);
	end = clock();
	cout << "Liczba iteracji: " << iteracje << " czas: " << double(end - begin) / CLOCKS_PER_SEC * 1000 << endl;
	//cout << (double)(end - begin) << endl;



	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);
	begin = clock();
	iteracje = solver.algorytmGaussaSeidela(macierz, wektor, x, N, BLAD);
	end = clock();
	cout << "Liczba iteracji: " << iteracje << " czas: " << double(end - begin) << endl;
	//cout << (double)(end - begin) << endl;
	delete[] x;
	delete[] wektor;
	delete[] macierz;
	system("PAUSE");
	return 0;
}
