// LinearEquationSolver.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "LinearEquationSolver.h"
#include <cmath>


// This is an example of an exported variable
LINEAREQUATIONSOLVER_API int nLinearEquationSolver=0;

// This is an example of an exported function.
LINEAREQUATIONSOLVER_API int fnLinearEquationSolver(void)
{
    return 42;
}

// This is the constructor of a class that has been exported.
// see LinearEquationSolver.h for the class definition
CLinearEquationSolver::CLinearEquationSolver()
{
    return;
}

double* CLinearEquationSolver::algorytmGaussa(double* macierz, double* wektor, double* x, int n)
{
	for (int h = 1; h < n; h++)
	{
		for (int i = h; i < n; i++)
		{
			double m = macierz[i*n + h - 1] / macierz[(h - 1)*n + h - 1];
			macierz[i*n + h - 1] = 0;
			for (int j = h; j < n; j++)
			{
				macierz[i*n + j] = macierz[i*n + j] - m*macierz[(h - 1)*n + j];
			}
			wektor[i] = wektor[i] - m*wektor[h - 1];
		}
	}

	x[n - 1] = wektor[n - 1] / macierz[(n - 1)*n + n - 1];

	for (int i = n - 2; i >= 0; i--) {
		double suma = 0;
		for (int j = i + 1; j < n; j++)
		{
			suma += macierz[i*n + j] * x[j];
		}
		x[i] = (wektor[i] - suma) / macierz[i*n + i];
	}
	return x;
}

double CLinearEquationSolver::stopienBledu(double* x, double* macierz, double* wektor, int n)
{
	double* wynik = new double[n];
	for (int i = 0; i < n; i++)
	{
		double suma = 0;
		for (int j = 0; j < n; j++)
		{
			suma += x[j] * macierz[i*n + j];
		}
		wynik[i] = suma;
	}
	double suma = 0;
	for (int i = 0; i < n; i++)
	{
		suma += pow(wynik[i] - wektor[i], 2);
	}
	return sqrt(suma);
}


int CLinearEquationSolver::algorytmJacobiego(double *macierz, double* wektor, double* x, int n, double blad)
{
	int iteracje = 0;
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}
	double *x_nowe = new double[n];

	while (stopienBledu(x, macierz, wektor, n) > blad)
	{
		for (int i = 0; i < n; i++)
		{
			double suma = 0;
			for (int j = 0; j < n; j++)
			{
				if (j == i)
				{
					continue;
				}
				suma -= macierz[i*n + j] * x[j];
			}
			x_nowe[i] = (suma + wektor[i]) / macierz[i*n + i];
		}
		for (int i = 0; i < n; i++)
		{
			x[i] = x_nowe[i];
		}
		iteracje++;
	}

	return iteracje;
}


int CLinearEquationSolver::algorytmGaussaSeidela(double *macierz, double *wektor, double *x, int n, double blad)
{
	int iteracje = 0;
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
	}
	double *x_nowe = new double[n];

	while (stopienBledu(x, macierz, wektor, n) > blad)
	{
		for (int i = 0; i < n; i++)
		{
			double suma = 0;
			for (int j = 0; j < i; j++)
			{
				suma -= macierz[i*n + j] * x_nowe[j];
			}
			double suma2 = 0;
			for (int j = i + 1; j < n; j++)
			{
				suma2 -= macierz[i*n + j] * x[j];
			}
			x_nowe[i] = (suma + suma2 + wektor[i]) / macierz[i*n + i];
		}
		for (int i = 0; i < n; i++)
		{
			x[i] = x_nowe[i];
		}
		iteracje++;
	}

	return iteracje;
}