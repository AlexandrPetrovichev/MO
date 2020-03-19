#include <vector>
#define test3
using namespace std;

typedef vector<vector<double>> Matrix;


double F1(double* x);
double GradF1(double* x);
double GradF2(double* x);
double GradFQuad(double* x);
double FindMin(double* s, double* p);
void PrintSolution(double* x, double val, int numIter);
double FletcherRievesMethod();
double inner_prod(double* x, double* y);

//for Pirson

double FindMin2(double* s, double* p);
double F(vector<double> &xk);
vector<double> dF(vector< double> xk);
vector<double>Sum(vector<double>& a, vector<double>& b);
vector<double> dG(vector<double>& xprev, vector<double>& xnext);
vector<double> Pirson(Matrix h0, vector<double> x0, double eps);
vector<double> MultMV(Matrix& hk, vector<double>& x);
vector<double> MultVM(Matrix& hk, vector<double>& x);
vector<double> MultMV(Matrix& hk, vector<double>& x);
Matrix MultvvT(vector<double>& a, vector<double>& b);
double MultvTv(vector<double>& a, vector<double>& b);
