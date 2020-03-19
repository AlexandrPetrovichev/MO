#include "sol.h"
#include <cmath>

typedef vector<vector<double>> Matrix;

vector<double>Sum(vector<double>& a, vector<double>& b)
{
	vector<double> result;
	for (size_t i = 0; i < a.size(); i++)
	{
		result.push_back(a[i] + b[i]);
	}
	return result;
}


vector<double> dG(vector<double>& xprev, vector<double>& xnext)
{
	int n = xprev.size();
	vector<double> Gprev = dF(xprev);
	vector<double>Gnext = dF(xnext);
	for (size_t i = 0; i < n; i++)
	{
		Gnext[i] = Gnext[i] - Gprev[i];
	}
	return Gnext;
}

vector<double> MultMV(Matrix& hk, vector<double>& x)
{
	vector<double> result;
	for (size_t i = 0; i < hk.size(); i++)
	{
		double sum = 0;
		for (size_t j = 0; j < hk[i].size(); j++)
		{
			sum += hk[i][j] * x[j];
		}
		result.push_back(sum);
	}
	return result;
}

vector<double> MultVM(Matrix& hk, vector<double>& x)
{
	vector<double> result;
	for (size_t i = 0; i < hk.size(); i++)
	{
		double sum = 0;
		for (size_t j = 0; j < hk[i].size(); j++)
		{
			sum += hk[j][i] * x[j];
		}
		result.push_back(sum);
	}
	return result;
}


Matrix MultvvT(vector<double>& a, vector<double>& b)
{
	int n = a.size();
	Matrix result(n);
	for (size_t i = 0; i < n; i++)
	{
		result[i].resize(n);
		for (size_t j = 0; j < n; j++)
		{
			result[i][j] = a[i] * b[j];
		}
	}
	return result;
}

double MultvTv(vector<double>& a, vector<double>& b)
{
	double result = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

vector<double> ResV(vector<double>& a, vector<double>& b)
{
	vector<double> res;
	for (size_t i = 0; i < a.size(); i++)
	{
		res.push_back(a[i] - b[i]);
	}
	return res;
}

double FindMin2(double* s, double* p)
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

vector<double> FindP(vector<double>& xk, vector<vector<double>>& h)
{
	vector<double> df = dF(xk);
	vector<double>tmp = MultMV(h,df);
	for (size_t i = 0; i < tmp.size(); i++)
	{
		tmp[i] *= -1;
	}
	return tmp;
}

Matrix Findh(Matrix& h,vector<double>& dx, vector<double> & dg)
{
	int n = h.size();
	Matrix dh;
	vector<double> tmp1;
	vector<double> tmp2;
	vector<double> tmp3;
	double tmp4;
	tmp1 = MultMV(h, dg);
	tmp1 = ResV(dx, tmp1);
	tmp2 = MultMV(h, dg);

	tmp3 = MultVM(h, dg);
	tmp4 = MultvTv(dg, tmp3);
	dh = MultvvT(tmp1,tmp2);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			dh[i][j] = h[i][j] + dh[i][j] / tmp4;
		}
	}
	return dh;
}

double Norm(vector<double>& x)
{
	double sum = 0;
	for (size_t i = 0; i < x.size(); i++)
	{
		sum += x[i] * x[i];
	}
	return sqrt(sum);
}

vector<double> Pirson(Matrix h0,
	vector<double> x0,double eps)
{
	int n = x0.size();
	vector<double> df,xpr,xk = x0,p,dg(n);
	double Lk;
	Matrix h = h0;
	int iter = 0;
	vector<double> dx(n);
	df = dF(xk);
	do
	{	
		iter++;
		p = FindP(xk, h);
		Lk = FindMin(xk.data(), p.data());
		xpr = xk;
		for (size_t i = 0; i < n; i++)
		{
			dx[i] = Lk * p[i];
			xk[i] += dx[i];
		}
		
		dg = dG(xpr, xk);
		h = Findh(h, dx, dg);
		df = dF(xk);
	} while (Norm(df)>eps);
	
	PrintSolution(xk.data(), F1(xk.data()), iter);
	return xk;
}
