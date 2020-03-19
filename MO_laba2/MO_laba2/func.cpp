#include <cmath>
// ���� �������
double F1(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = x[0] * x[0] + x[1] * x[1] - 5;
    //res = 100 * (x[0] - x[1]) * (x[0] - x[1]) + (x[0] - 1) * (x[0] - 1);
    res = 2 * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2) +
        3 * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    return -res;
}
// ����������� �� ������ ����������
double GradF1(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[0];
    //res = 200 * (x[0] - x[1]) + 2 * (x[0] - 1);
    res = -(x[0] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[0] - 1) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    return -res;
}
// ����������� �� ������ ����������
double GradF2(double* x)
{
    double a = 0.1, eps = 0.00001, res = 0;
    //res = 2 * x[1];
    res = -(x[1] - 1) * exp(-(x[0] - 1) / 2 * (x[0] - 1) / 2 - (x[1] - 1) / 2 * (x[1] - 1) / 2)
        - 2 / 3 * (x[1] - 3) * exp(-(x[0] - 1) / 3 * (x[0] - 1) / 3 - (x[1] - 3) / 3 * (x[1] - 3) / 3);
    //res = -200 * (x[0] - x[1]);
    return -res;
}