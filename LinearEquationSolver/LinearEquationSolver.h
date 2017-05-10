// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the LINEAREQUATIONSOLVER_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// LINEAREQUATIONSOLVER_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef LINEAREQUATIONSOLVER_EXPORTS
#define LINEAREQUATIONSOLVER_API __declspec(dllexport)
#else
#define LINEAREQUATIONSOLVER_API __declspec(dllimport)
#endif

// This class is exported from the LinearEquationSolver.dll
class LINEAREQUATIONSOLVER_API CLinearEquationSolver {
public:
	CLinearEquationSolver(void);
	double* algorytmGaussa(double* macierz, double* wektor, double* x, int n);
	double stopienBledu(double* x, double* macierz, double* wektor, int n);
	int algorytmJacobiego(double *macierz, double* wektor, double* x, int n, double blad);
	int algorytmGaussaSeidela(double *macierz, double *wektor, double *x, int n, double blad);
	// TODO: add your methods here.
};

extern LINEAREQUATIONSOLVER_API int nLinearEquationSolver;

LINEAREQUATIONSOLVER_API int fnLinearEquationSolver(void);
