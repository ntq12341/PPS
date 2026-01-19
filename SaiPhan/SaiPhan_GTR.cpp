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

double r(double x){
	return 1 + cos(x)*cos(x);
}

matrix SP_VTR(double x0, double xn, double h){
	int n = (xn - x0)/h; 
	matrix M(n-1, n-1);
	matrix D(n-1, n-1);

	double a, b, c; 
	for(int i = 1;i < n;i++){
		a = p(x0 + h*(i - 0.5)); 
		b = p(x0 + h*(i + 0.5)) + p(x0 + h*(i - 0.5)) + h*h*q(x0 + h*i);
		c = p(x0 + h*(i + 0.5)); 
		if(i == 1){
			M.mat[i-1][i-1] = -b;
			M.mat[i-1][i] = c;
		}else if(i == n-1){
			M.mat[i-1][i-2] = a;
			M.mat[i-1][i-1] = -b;
		}else{
			M.mat[i-1][i-2] = a;
			M.mat[i-1][i-1] = -b;
			M.mat[i-1][i] = c; 
		}
		D.mat[i-1][i-1] = -h*h*r(x0 + i*h);
	}

	matrix D1 = D.inverse();
	matrix U = D1*M;
	// M.out();
	// D1.out();
	// U.out();
	return U; 
}

int main(){
	double x0 = 0, xn = 3, h = 0.1;
	matrix U = SP_VTR(x0, xn, h);
	U.out();
}