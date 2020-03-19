#include <iostream>
#include <cstdlib>  
#include <cmath>
#include "sol.h"
using namespace std;

double FletcherRievesMethod()
{
    double eps = 1e-3;
    const double EPS = 1e-5;
    //Начальное приближение
    double x[2] = { 1, 2 };
    double curVal = F1(x);
    double prevVal = curVal;
    double p[2];
    p[0] = -GradF1(x);
    p[1] = -GradF2(x);

    double gradSquare = inner_prod(p, p);

    int numIter = 0;
    do
    {
        numIter++;
        double lambda, beta, newGradSquare;
        double newGrad[2];


        //Ищем минимум F1(x + lambda * p) с помощью метода одномерной оптимизации
        lambda = FindMin(x, p);
        x[0] = x[0] + lambda * p[0];
        x[1] = x[1] + lambda * p[1];
        newGrad[0] = -GradF1(x);
        newGrad[1] = -GradF2(x);

        newGradSquare = inner_prod(newGrad, newGrad);

        //if (numIter % (10) == 0) beta = 0; //Обновление
        //else
            beta = newGradSquare / gradSquare; //Используем метод Флетчера - Ривса
        p[0] = newGrad[0] + beta * p[0];
        p[1] = newGrad[1] + beta * p[1];

        prevVal = curVal;
        curVal = F1(x);

        gradSquare = newGradSquare;
    } while (gradSquare > eps);

    cout << "Fletcher-Rieves Method: " << endl;
    PrintSolution(x, F1(x), numIter);
    return 0;

}