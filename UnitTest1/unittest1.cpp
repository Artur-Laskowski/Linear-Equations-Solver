#include "stdafx.h"
#include "CppUnitTest.h"
#include "../LinearEquationSolver/LinearEquationSolver.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			Assert::AreEqual(1, 1);
		}

		TEST_METHOD(gaussianTest)
		{
			//Arrange
			double x[4];
			int n = 4;

			double wektor[] = { 13.15, 49.84, -14.08, -46.51 };

			double macierz[] = {
				1.2, 2.6, -0.1, 1.5,
				4.5, 9.8, -0.4, 5.7,
				0.1, -0.1, -0.3, -3.5,
				4.5, -5.2, 4.2, -3.4,
			};

			//Act
			CLinearEquationSolver linearSolver;
			linearSolver.algorytmGaussa(macierz, wektor, x, n);

			//Assert
			double wyniki[4] = {-1.3, 3.2, -2.4, 4.1};
			for (int i = 0; i < n; i++)
			{
				Assert::AreEqual(wyniki[i], x[i], 10e-4);
			}
		}

		TEST_METHOD(stopienBleduTest)
		{
			//Arrange
			double x[4];
			int n = 4;

			double wektor[] = { 13.15, 49.84, -14.08, -46.51 };

			double macierz[] = {
				1.2, 2.6, -0.1, 1.5,
				4.5, 9.8, -0.4, 5.7,
				0.1, -0.1, -0.3, -3.5,
				4.5, -5.2, 4.2, -3.4,
			};

			double wyniki[4] = { -1.3, 3.2, -2.4, 4.1 };

			//Act
			CLinearEquationSolver linearSolver;
			double blad = linearSolver.stopienBledu(wyniki, macierz, wektor, n);

			//Assert
			Assert::AreEqual(0.0, blad, 10e-9);
		}

		TEST_METHOD(jacobianTest)
		{
			//Arrange
			double x[4];
			int n = 4;

			double wektor[] = { 6, 25, -11, 15 };

			double macierz[] = {
				10, -1, 2, 0,
				-1, 11, -1, 3,
				2, -1, 10, -1,
				0, 3, -1, 8,
			};

			//Act
			CLinearEquationSolver linearSolver;
			int iteracje = linearSolver.algorytmJacobiego(macierz, wektor, x, n, 10e-9);

			//Assert

			double wyniki[4] = { 1, 2, -1, 1 };
			for (int i = 0; i < n; i++)
			{
				Assert::AreEqual(wyniki[i], x[i], 10e-4);
			}
		}

		TEST_METHOD(gaussSeidelTest)
		{
			//Arrange
			double x[4];
			int n = 4;

			double wektor[] = { 6, 25, -11, 15 };

			double macierz[] = {
				10, -1, 2, 0,
				-1, 11, -1, 3,
				2, -1, 10, -1,
				0, 3, -1, 8,
			};

			//Act
			CLinearEquationSolver linearSolver;
			linearSolver.algorytmGaussaSeidela(macierz, wektor, x, n, 10e-9);

			//Assert
			double wyniki[4] = { 1, 2, -1, 1 };
			for (int i = 0; i < n; i++)
			{
				Assert::AreEqual(wyniki[i], x[i], 10e-5);
			}
		}
	};
}