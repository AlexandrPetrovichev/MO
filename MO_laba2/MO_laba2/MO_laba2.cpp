#include <iostream>
#include <cstdlib>  
#include <cmath>
#include "sol.h"
using namespace std;

// Вывод результата на экран
void PrintSolution(double* x, double val, int numIter)
{
    cout << "-----------------------" << endl;
    cout << "Number of iterations: " << numIter << endl;
    cout << "Computed solution: " << endl << "[" << x[0] << "," << x[1] << "]" << endl;
    cout << "Function value: " << val << endl << endl;
}

// Вычисление скалярного произведения
double inner_prod(double* x, double* y) { return x[1] * y[1] + x[0] * y[0] - 5; }

int main()
{
    cout.precision(9);
    FletcherRievesMethod();
    return 0;
}


/*

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <conio.h>


using namespace std;
double f(double *x)
{
    return x[0] * x[0] + x[1] * x[1];
    //return 100 * (x[0] - x[1]) * (x[0] - x[1]) + (x[0] - 1) * (x[0] - 1);
}
double* Nabla(double* x)
{
	double nabla[2];
	//nabla[0] = (200 * (x[0] - x[1]) + 2 * (x[0] - 1));
	//nabla[1] = (-200 * (x[0] - x[1]));
    nabla[0] = 2 * x[0];
    nabla[1] = 2 * x[1];
	return nabla;
}

double Golden(double* x, double* S)
{
    double a = -10, b = 10, eps = 0.00001;
    double r = (3 - sqrt(5)) / 2;
    double f1, f2;
    double x1 = a + r * (b - a);
    double x2 = b - r * (b - a);
    double y1[2], y2[2];
    y1[0] = (x[0] + x1 * S[0]);
    y1[1] = (x[1] + x1 * S[1]);
    y2[0] = (x[0] + x2 * S[0]);
    y2[1] = (x[1] + x2 * S[1]);
    f1 = f(y1);    
    f2 = f(y2);
    if (f1 < f2)
    {
        b = x2;
        x2 = x1;
        y2[0] = y1[0];
        y2[1] = y1[1];
        x1 = a + r * (b - a);
        y1[0] = x[0] + x1 * S[0];
        y1[1] = x[1] + x1 * S[1];
        
        f1 = f(y1);
    }
    else
    {
        a = x1;
        x1 = x2;
        y1[0] = y2[0];
        y1[1] = y2[1];
        x2 = b - r * (b - a);
        y2[0] = x[0] + x2 * S[0];
        y2[1] = x[1] + x2 * S[1];
        
        f2 = f(y2);
    }
    while (abs(a - b) > eps)
    {
        if (f1 < f2)
        {
            b = x2;
            x2 = x1;
            y2[0] = y1[0];
            y2[1] = y1[1];
            x1 = a + r * (b - a);
            y1[0] = x[0] + x1 * S[0];
            y1[1] = x[1] + x1 * S[1];
            
            f1 = f(y1);
        }
        else
        {
            a = x1;
            x1 = x2;
            y1[0] = y2[0];
            y1[1] = y2[1];
            x2 = b - r * (b - a);
            y2[0] = x[0] + x2 * S[0];
            y2[1] = x[1] + x2 * S[1];
           
            f2 = f(y2);
        }
    }
    return (a + b) / 2;
}





double Norma(double* x)
{
	double norma = 0;
	for (int i = 0; i < 2; i++)
		norma = x[i] * x[i];
	norma = sqrt(norma);
	return norma;
}
void MSG(double* x_k, double eps, double lambda)
{
	
	double S[2], x_k1[2], SS[2];
    S[0] = -Nabla(x_k)[0];
    S[1] = -Nabla(x_k)[1];
    for (int i = 0; i < 2; i++)
        SS[i] = S[i];
	int iter = 1;
	double lambda_k;
    while (Norma(SS) >= eps)
    {
        for (int i = 0; i < 2; i++)
            SS[i] = S[i];
        for (int i = 0; i < 2; i++)
            x_k1[i] = x_k[i];
        double arg[2];
        lambda_k = Golden(x_k, S);
        for (int i = 0; i < 2; i++)
            x_k[i] += lambda_k * S[i];
        for (int i = 0; i < 2; i++)
        {
            double w = (Norma(Nabla(x_k1)) * Norma(Nabla(x_k1))) / (Norma(Nabla(x_k)) * Norma(Nabla(x_k)));
            S[i] = S[i] * w - Nabla(x_k)[i] * x_k[i];
        } 
        printf_s("%f\t%f\n", x_k[0], x_k[1]);
	}
   
}
void main()
{
	double x[2];
    x[0] = 0.5;
    x[1] = 0.5;
	MSG(x, 1e-7, 1);

}

*/