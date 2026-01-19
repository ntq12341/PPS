#include "..\Gauss-Jordan\Gauss-Jordan.cpp"
using namespace std;
#define PI 3.14159265359
#define E 2.718281828


double p(double x){
	return x*x + 1; 
}

double q(double x){
	return 4*pow(x,4) - 4*x*x + 5;
}

double f(double x){
	double e = E;
	return 1.25*x*pow(e,-x*x); 
}

vector<double> SaiPhan1(double x0, double xn, double h, double a, double b){
	int n = (xn - x0)/h; 
	matrix M(n+1, n+1);
	matrix B(n+1,1);

	M.mat[0][0] = 1;
	M.mat[n][n] = 1;

	B.mat[0][0] = a;
	B.mat[n][0] = b;
	for(int i = 1;i < n;i++){
		M.mat[i][i-1] = p(x0 + (i - 0.5) * h);
		M.mat[i][i] = -(p(x0 + (i + 0.5) * h) + p(x0 + (i - 0.5) * h) + h*h*q(x0 + i*h)) ;
		M.mat[i][i+1] = p(x0 + (i + 0.5) * h);

		B.mat[i][0] = -h*h*f(x0 + i*h);
	}

	matrix U = gauss_jordan(M,B);

	vector<double> u = U.t().mat[0];

	return u;
}

vector<double> SaiPhan2(double x0, double xn, double h, double a, double b){
	int n = (xn - x0)/h; 
	matrix M(n+1, n+1);
	matrix B(n+1,1);

	M.mat[0][0] = -(p(x0 + 0.5*h) + h*h*q(x0)*0.5);
	M.mat[0][1] = p(x0 + 0.5*h); 
	M.mat[n][n-1] = -p(x0 + (n-0.5)*h);
	M.mat[n][n] = p(x0 + (n-0.5)*h) + h*h*q(x0 + n*h)*0.5;

	B.mat[0][0] = -h*h*f(x0)*0.5 - a*h;
	B.mat[n][0] = h*h*f(x0 + n*h)*0.5 - b*h;
	for(int i = 1;i < n;i++){
		M.mat[i][i-1] = p(x0 + (i - 0.5) * h);
		M.mat[i][i] = -(p(x0 + (i + 0.5) * h) + p(x0 + (i - 0.5) * h) + h*h*q(x0 + i*h)) ;
		M.mat[i][i+1] = p(x0 + (i + 0.5) * h);

		B.mat[i][0] = -h*h*f(x0 + i*h);
	}

	matrix U = gauss_jordan(M,B);

	vector<double> u = U.t().mat[0];

	return u;
}

vector<double> SaiPhan3(double x0, double xn, double h, double a, double b, double o1, double o2){
	int n = (xn - x0)/h; 
	matrix M(n+1, n+1);
	matrix B(n+1,1);

	M.mat[0][0] = -(p(x0 + 0.5*h) + h*h*q(x0)*0.5 + o1);
	M.mat[0][1] = p(x0 + 0.5*h); 
	M.mat[n][n-1] = -p(x0 + (n-0.5)*h);
	M.mat[n][n] = p(x0 + (n-0.5)*h) + h*h*q(x0 + n*h)*0.5 - o2;

	B.mat[0][0] = -h*h*f(x0)*0.5 - a*h;
	B.mat[n][0] = h*h*f(x0 + n*h)*0.5 - b*h;
	for(int i = 1;i < n;i++){
		M.mat[i][i-1] = p(x0 + (i - 0.5) * h);
		M.mat[i][i] = -(p(x0 + (i + 0.5) * h) + p(x0 + (i - 0.5) * h) + h*h*q(x0 + i*h)) ;
		M.mat[i][i+1] = p(x0 + (i + 0.5) * h);

		B.mat[i][0] = -h*h*f(x0 + i*h);
	}

	matrix U = gauss_jordan(M,B);

	vector<double> u = U.t().mat[0];

	return u;
}