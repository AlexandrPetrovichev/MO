#include <iostream>
#include <cstdlib>  
#include <cmath>
#include "sol.h"

//Одномерная оптимизация - применим метод золотого сечения
double FindMin(double* s, double* p)
{
    const double eps = 1e-8;
    const double r = 1.618;
    double a = 0;
    double b = 10;
    double x0, x1, xf1, xf2;
    x0 = b - (b - a) / r; 
    x1 = a + (b - a) / r; 
    double t1[2], t2[2];
    do {
        t1[0] = s[0] + x0 * p[0];
        t1[1] = s[1] + x0 * p[1];
        t2[0] = s[0] + x1 * p[0];
        t2[1] = s[1] + x1 * p[1];
        xf1 = F1(t1);
        xf2 = F1(t2); 
        if (xf1 >= xf2)
        {
            a = x0;
            x0 = x1;
            t2[0] = s[0] + x1 * p[0];
            t2[1] = s[1] + x1 * p[1];
            xf1 = F1(t2);
            t2[0] = s[0] + x1 * p[0];
            t2[1] = s[1] + x1 * p[1];
            x1 = a + (b - a) / r;
            xf2 = F1(t2);
        }
        else
        {
            b = x1;
            x1 = x0;
            xf2 = xf1;
            x0 = b - (b - a) / r;
            t1[0] = s[0] + x0 * p[0];
            t1[1] = s[1] + x0 * p[1];
            xf1 = F1(t1);
        }
    } while (fabs(b - a) > eps);
    return (a + b) / 2;
}