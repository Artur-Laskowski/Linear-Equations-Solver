#include <iostream>
#include <math.h>
#include <ctime>

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
}

void print(double* macierz)
{
	for (int i = 0; i < N; i++)
	{
		cout << macierz[i] << ", ";
	}
	cout << endl;
}

double stopienBledu(double* x, double* macierz, double* wektor)
{
	double* wynik = new double[N];
	for (int i = 0; i < N; i++)
	{
		double suma = 0;
		for (int j = 0; j < N; j++)
		{
			suma += x[j] * macierz[i*N + j];
		}
		wynik[i] = suma;
	}
	double suma = 0;
	for (int i = 0; i < N; i++)
	{
		suma += pow(wynik[i] - wektor[i], 2);
	}
	return sqrt(suma);
}

double* algorytmGaussa(double* macierz, double* wektor, double* x)
{
	for (int h = 1; h < N; h++)
	{
		for (int i = h; i < N; i++)
		{
			double m = macierz[i*N + h - 1] / macierz[(h - 1)*N + h - 1];
			macierz[i*N + h - 1] = 0;
			for (int j = h; j < N; j++)
			{
				macierz[i*N + j] = macierz[i*N + j] - m*macierz[(h-1)*N + j];
			}
			wektor[i] = wektor[i] - m*wektor[h-1];
		}
	}

	x[N-1] = wektor[N - 1] / macierz[(N - 1)*N + N - 1];

	for (int i = N-2; i >= 0; i--) {
		double suma = 0;
		for (int j = i+1; j < N; j++)
		{
			suma += macierz[i*N + j] * x[j];
		}
		x[i] = (wektor[i] - suma) / macierz[i*N + i];
	}
	return x;
}


int algorytmJacobiego(double *macierz, double* wektor, double* x)
{
	int iteracje = 0;
	for (int i = 0; i < N; i++)
	{
		x[i] = 0;
	}
	double *x_nowe = new double[N];

	while (stopienBledu(x, macierz, wektor) > BLAD)
	{
		for (int i = 0; i < N; i++)
		{
			double suma = 0;
			for (int j = 0; j < N; j++)
			{
				if (j == i)
				{
					continue;
				}
				suma -= macierz[i*N + j] * x[j];
			}
			x_nowe[i] = (suma + wektor[i]) / macierz[i*N + i];
		}
		for (int i = 0; i < N; i++)
		{
			x[i] = x_nowe[i];
		}
		iteracje++;
	}

	return iteracje;
}


int algorytmGaussaSeidela(double *macierz, double *wektor, double *x)
{
	int iteracje = 0;
	for (int i = 0; i < N; i++)
	{
		x[i] = 0;
	}
	double *x_nowe = new double[N];

	while (stopienBledu(x, macierz, wektor) > BLAD)
	{
		for (int i = 0; i < N; i++)
		{
			double suma = 0;
			for (int j = 0; j < i; j++)
			{
				suma -= macierz[i*N + j] * x_nowe[j];
			}
			double suma2 = 0;
			for (int j = i + 1; j < N; j++)
			{
				suma2 -= macierz[i*N + j] * x[j];
			}
			x_nowe[i] = (suma + suma2 + wektor[i]) / macierz[i*N + i];
		}
		for (int i = 0; i < N; i++)
		{
			x[i] = x_nowe[i];
		}
		iteracje++;
	}

	return iteracje;
}


void wypelnijWektor(double* wektor)
{
	for (int i = 0; i < N; i++)
	{
		wektor[i] = sin(static_cast<float>(i + 1)*(F + 1) / 50);
	}
}

int main()
{
	auto macierz = new double[N*N];
	auto x = new double[N];
	auto wektor = new double[N];

#ifdef DEBUG
	wektor[0] = 13.15;
	wektor[1] = 49.84;
	wektor[2] = -14.08;
	wektor[3] = -46.51;

	macierz[0] = 1.2;
	macierz[1] = 2.6;
	macierz[2] = -0.1;
	macierz[3] = 1.5;
	macierz[4] = 4.5;
	macierz[5] = 9.8;
	macierz[6] = -0.4;
	macierz[7] = 5.7;
	macierz[8] = 0.1;
	macierz[9] = -0.1;
	macierz[10] = -0.3;
	macierz[11] = -3.5;
	macierz[12] = 4.5;
	macierz[13] = -5.2;
	macierz[14] = 4.2;
	macierz[15] = -3.4;
#else
	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);
#endif


	clock_t begin = clock();
	algorytmGaussa(macierz, wektor, x);
	cout << "Norma z residuum: " << stopienBledu(x, macierz, wektor) << endl;
	clock_t end = clock();
	cout << "Czas: " << double(end - begin) / CLOCKS_PER_SEC * 1000 << endl;
	//cout << (double)(end - begin) << endl;
	//print(x);

	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);
	begin = clock();
	int iteracje = algorytmJacobiego(macierz, wektor, x);
	end = clock();
	cout << "Liczba iteracji: " << iteracje << " czas: " << double(end - begin) / CLOCKS_PER_SEC * 1000 << endl;
	//cout << (double)(end - begin) << endl;



	wypelnijMacierz(macierz);
	wypelnijWektor(wektor);
	begin = clock();
	iteracje = algorytmGaussaSeidela(macierz, wektor, x);
	end = clock();
	cout << "Liczba iteracji: " << iteracje << " czas: " << double(end - begin) << endl;
	//cout << (double)(end - begin) << endl;
	delete[] x;
	delete[] wektor;
	delete[] macierz;
	system("PAUSE");
	return 0;
}