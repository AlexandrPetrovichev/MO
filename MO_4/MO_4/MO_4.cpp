#include <iostream>
#include <math.h>
#include <random>
#include <omp.h>
#include <fstream>
using namespace std;
double eps = 1e-7;
int N = 1;
double fc1 = 0, fc2 = 0, fc3 = 0;
double m = 100;
double r = 1;
double A[2] = { -10, -10 };
double B[2] = { 10, 10 };
double Norm(double* x, double* y) { return sqrt(x[0] * y[0] + x[1] * y[1]); }
double Norm(pair<double, double> x, pair<double, double> y) { return sqrt(x.first * y.first + x.second * y.second); }
pair<double, double> Rand()
{
	pair<double, double> x;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-10, 10);
	x.first = dis(gen);
	x.second = dis(gen);
	return x;
}
double g(pair<double, double> x){ return -x.first - 10; }
double h(pair<double, double> x) { return x.first - 10; }
double p(pair<double, double> x) { return -x.second - 10; }
double q(pair<double, double> x) { return x.second - 10; }
double G(pair<double, double> x) { return pow(0.5 * (g(x) + abs(g(x))), 2); }
double H(pair<double, double> x) { return pow(0.5 * (h(x) + abs(h(x))), 2); }
double P(pair<double, double> x) { return pow(0.5 * (p(x) + abs(p(x))), 2); }
double Q(pair<double, double> x) { return pow(0.5 * (q(x) + abs(q(x))), 2); }
double r_correct(double r) { return r * 2; }
double f(pair<double, double> z)
{
	double f = 0;
	double C[6] = { 2,3,8,3,2,8 };
	double a[6] = { 3,-5,0,3,-4,6 };
	double b[6] = { -4,-6,-1,7,0,5 };
	double x = z.first, y = z.second;
	for (int i = 0; i < 6; i++)
	{
		f += C[i] / (1 + (x - a[i]) * (x - a[i]) + (y - b[i]) * (y - b[i]));
	}
	f -= r * (G(z) + H(z) + P(z) + Q(z));
	return -f;
}
double ff(pair<double, double> z)
{
	double f = 0;
	double C[6] = { 2,3,8,3,2,8 };
	double a[6] = { 3,-5,0,3,-4,6 };
	double b[6] = { -4,-6,-1,7,0,5 };
	double x = z.first, y = z.second;
	for (int i = 0; i < 6; i++)
	{
		f += C[i] / (1 + (x - a[i]) * (x - a[i]) + (y - b[i]) * (y - b[i]));
	}
	fc1++;
	fc2++;
	fc3++;
	return -f;
}

pair<double, double> Easy_random_search()
{
	pair<double, double> x, x1;
	int k = 0;
	x = Rand();
	for (int i = 0; i < N; i++)
	{
		x1 = Rand();
		if (ff(x1) < ff(x))
			x = x1;
	}
	//cout << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	return x;
}
double Golden(pair<double, double> s, double* p, int& fc)
{
	const double r = (1 + sqrt(5)) / 2;
	double a = -10;
	double b = 10;
	double x0, x1, xf1, xf2;
	x0 = b - (b - a) / r;
	x1 = a + (b - a) / r;
	pair<double, double> t1, t2;
	t1.first = s.first + x0 * p[0];
	t1.second = s.second + x0 * p[1];
	t2.first = s.first + x1 * p[0];
	t2.second = s.second + x1 * p[1];
	xf1 = f(t1);
	xf2 = f(t2);
	fc += 2;
	do {

		if (xf1 >= xf2)
		{
			a = x0;
			x0 = x1;
			xf1 = xf2;
			x1 = a + (b - a) / r;
			t2.first = s.first + x1 * p[0];
			t2.second = s.second + x1 * p[1];
			xf2 = f(t2);
			fc++;
		}
		else
		{
			b = x1;
			x1 = x0;
			xf2 = xf1;
			x0 = b - (b - a) / r;
			t1.first = s.first + x0 * p[0];
			t1.second = s.second + x0 * p[1];
			xf1 = f(t1);
			fc++;
		}
	} while (abs(b - a) > eps);
	return (a + b) / 2;
}
pair<double, double> Roz(pair<double, double> x)
{
	int fc = 0, i = 0;
	pair<double, double> x_old;
	double a1[2], a2[2];
	double s1[2] = { 1 , 0 }, s2[2] = { 0, 1 };
	double lambda1;
	double lambda2;
	do {

		x_old.first = x.first;
		x_old.second = x.second;
		lambda1 = Golden(x, s1, fc);
		x.first = x.first + s1[0] * lambda1;
		x.second = x.second + s1[1] * lambda1;
		lambda2 = Golden(x, s2, fc);
		x.first = x.first + s2[0] * lambda2;
		x.second = x.second + s2[1] * lambda2;

		a1[0] = s1[0] * lambda1 + s2[0] * lambda2;
		a1[1] = s1[1] * lambda1 + s2[1] * lambda2;
		a2[0] = s2[0] * lambda2;
		a2[1] = s2[1] * lambda2;

		s1[0] = a1[0] / Norm(a1, a1);
		s1[1] = a1[1] / Norm(a1, a1);

		s2[0] = (a2[0] * Norm(a1, a1) * Norm(a1, a1) - a1[0] * Norm(a2, a2) * Norm(a2, a2)) /
			(Norm(a1, a1) * Norm(a2, a2) * sqrt(Norm(a1, a1) * Norm(a1, a1) - Norm(a2, a2) * Norm(a2, a2)));
		s2[1] = (a2[1] * Norm(a1, a1) * Norm(a1, a1) - a1[1] * Norm(a2, a2) * Norm(a2, a2)) /
			(Norm(a1, a1) * Norm(a2, a2) * sqrt(Norm(a1, a1) * Norm(a1, a1) - Norm(a2, a2) * Norm(a2, a2)));
		i++;
		r = r_correct(r);
	} while (abs(f(x_old) - f(x)) > eps);
	return x;
}
pair<double, double> Random_vec()
{
	pair<double, double> x;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-0.1, 0.1);
	x.first = dis(gen);
	x.second = sqrt(1 - x.first * x.first);
	return x;
}
pair<double, double> Global_search_1()
{
	pair<double, double> x, x1;
	int flag = 0, q = 0;
	x = { 5, 5 };
	while (flag < m)
	{
		x1 = Rand();
		x1 = Roz(x1);
		if (ff(x1) < ff(x))
		{
			x = x1;
			q = 1;
		}
		else
		{
			flag++;
			q = 0;
		}
		if (q == 1)
			flag = 0;
	}
	//cout << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	return x;
}

pair<double, double> Global_search_2()
{
	pair<double, double> x, x1;
	int flag = 0, q = 0, i = 0;;
	x = Rand();
	x = Roz(x);
	do {
		flag = 0;
		do
		{
			x1 = Easy_random_search();
			flag++;
		} while (ff(x1) > ff(x) && flag < m);
		if (flag < m) {
			x1 = Roz(x1);
			x = Rand();
			x = Easy_random_search();
		}
		else
			q = 1;
		i++;
	} while (q != 1);
	return x;
}



pair<double, double>Best_prob(pair<double, double> x)
{
	pair<double, double> e, x1, x2;
	e = Random_vec();
	
	return x;
}
double Sign(double x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}
pair<double, double> Global_search_3()
{
	pair<double, double> x1 = {0,0 }, x2,e;
	int i = 0, f = m;
	do {
	x1 = Roz(x1);
	x2 = x1;
		do
		{
			e = Random_vec();
			x2.first += e.first;
			x2.second += e.second;
			i++;
		} while (ff(x2) <= ff(x1) && fabs(x2.first) < 10 && fabs(x2.second) < 10);
		if (fabs(x2.first) > 10 && fabs(x2.second) > 10)
		{
			x2.first -= e.first;
			x2.second -= e.second;
		}
		x2 = Roz(x2);
		if (ff(x2) < ff(x1))
			x1 = x2;
		else
			f--;
	} while (f > 0);
	return x2;
}
int main()
{
	pair<double, double> x;
	ofstream fout;
	fout.open("out.txt");
//	eps = 1;
//	for(int j = 0; j < 3; j++)
//	{
//		eps = 1 / pow(10, j);
//		double P = 0.9;
//		double V = (B[0] - A[0]) * (B[1] - A[1]);
//		double Veps = eps * eps;
//		double Peps = Veps / V;
//		double a = log(1 - P) / log(1 - Peps);
//		while (N < log(1 - P) / log(1 - Peps))
//			N += 100;
//		x = Easy_random_search();
//		fout << eps << '\t' << P << '\t' << N << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
//	}
	eps = 1e-2;
	for (int i = 0; i < 3; i++)
	{
		double P = 0.1 + 0.4 * i;
		double V = (B[0] - A[0]) * (B[1] - A[1]);
		double Veps = eps * eps;
		double Peps = Veps / V;
		double a = log(1 - P) / log(1 - Peps);
		while (N < log(1 - P) / log(1 - Peps))
			N += 100;
		x = Easy_random_search();;
		fout << eps << '\t' << P << '\t' << N << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;

	}
	//m = 1;
	//for (int i = 0; i < 3; i++)
	//{
	//	m *= 10;
	//	fc1 = 0;
	//	fc2 = 0;
	//	fc3 = 0;
	//	fout.precision(8);
	//	x = Global_search_1();
	//	fout << "1" << '\t';
	//	fout << m << '\t' << fc1 << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	//	x = Global_search_2();
	//	fout << "2" << '\t';
	//	fout << m << '\t' << fc2 << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	//	x = Global_search_3();
	//	fout << "3" << '\t';
	//	fout << m << '\t' << fc3 << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	//
	//}
	//m = 100;
	//
	//x = Global_search_1();
	//	fout << "1" << '\t';
	//	cout << m << '\t' << -x.first << '\t' << -x.second << '\t' << -ff(x) << endl;
	fout.close();

	return 0;
}


