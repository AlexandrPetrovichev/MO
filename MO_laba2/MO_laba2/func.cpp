#include <cmath>
#include "sol.h"
// Сама функция
using namespace std;

#ifdef test1
double F1(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = x[0] * x[0] + x[1] * x[1] - 5;
    //res = 100 * (x[0] - x[1]) * (x[0] - x[1]) + (x[0] - 1) * (x[0] - 1);
    res = 2 * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2) +
        3 * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    return -res;
}
// Производная по первой переменной
double GradF1(double* xk)
{
    double a = 0.1, eps = 0.00001, res = 0;
    double x = xk[0];
    double y = xk[1];
    res = 1.0 / 3 * exp(-1.0 / 4 * (x - 1) * (x - 1)) * (x - 1) * (-2.0 * exp(5.0 / 36 * (x - 1) * (x - 1)) -
        3 * exp(-1.0 / 4 * (y - 1) * (y - 1)));
    return -res;
}
// Производная по второй переменной
double GradF2(double* xk)
{
    double a = 0.1, eps = 0.00001, res = 0;
    double x=xk[0];
    double y = xk[1];
    res = -2.0/3*(y-3)*exp(-1.0/9*(x-1)* (x - 1)-1.0/9*(y-3)*(y-3))-
        (y-1)*exp(-1.0/4*(x-1)*(x-1)-1.0/4*(y-1)*(y-1));
    //res = -200 * (x[0] - x[1]);
    return -res;
}

vector<double> dF(vector< double> xk)
{
    double a = 0.1, eps = 0.00001, res = 0;
    vector < double> result;
    result.push_back(GradF1(xk.data()));
    result.push_back(GradF2(xk.data()));

    return result;
}

#endif

#ifdef test2
double F1(double* xk)
{
    double x = xk[0];
    double y = xk[1];
    return 100 * (y - x) * (y - x) + (1 - x) * (1 - x);
}
// Производная по первой переменной
double GradF1(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[0];
    //res = 200 * (x[0] - x[1]) + 2 * (x[0] - 1);
    res = -(x[0] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[0] - 1) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    return -res;
}
// Производная по второй переменной
double GradF2(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[1];
    res = -(x[1] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[1] - 3) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    //res = -200 * (x[0] - x[1]);
    return -res;
}

vector<double> dF(vector< double> xk)
{
    double x = xk[0];
    double y = xk[1];
    vector<double> result(2);

    return { -200 * (y - x) - 2 * (1 - x),200 * (y - x) };
}
#endif

#ifdef test3
double F1(double* xk)
{
    double x = xk[0];
    double y = xk[1];
    return 100 * (y - x * x) * (y - x * x) + (1 - x) * (1 - x);
}
// Производная по первой переменной
double GradF1(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[0];
    //res = 200 * (x[0] - x[1]) + 2 * (x[0] - 1);
    res = -(x[0] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[0] - 1) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    return -res;
}
// Производная по второй переменной
double GradF2(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[1];
    res = -(x[1] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[1] - 3) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    //res = -200 * (x[0] - x[1]);
    return -res;
}

vector<double> dF(vector< double> xk)
{
    double x = xk[0];
    double y = xk[1];
    vector<double> result(2);

    return { 2*(200*x*x*x-200*x*y+x-1),200 * (y - x*x) };
}
#endif