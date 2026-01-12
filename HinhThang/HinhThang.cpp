#include "../Inverse/Inverse.cpp"
#include <cmath>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
double f(double x);
double hinhthang_case1(double a,double b, double eps);
double hinhthang_case2(double a,double b, double eps);
double hinhthang_case3(vector<double> &X, vector<double> &Y, double &eps); 

//------------------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(){
    cout << hinhthang_case2(0,2,1e-6);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double f(double x){
	return 1/(x*x +1);
}

double hinhthang_case1(double a,double b, double eps){
	double M2 = 2;
	int n = floor((M2/12) * (pow(b-a, 3)/eps));
	double x0, x1, y0, y1;
	double ans = 0;
	for(int i = 0;i < n;i++){
		x0 = a + i*(b-a)/n;
		x1 = x0 + (b-a)/n;
		y0 = f(x0);
		y1 = f(x1);
		//cout << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1 << '\n';
		ans += y0 + y1;
	}
	ans *= (b-a)/(2*n);
	return ans ;
}

double hinhthang_case2(double a,double b, double eps){
	int n = 2;
	int p = 2;
	int q = n*p;
	double In, Ipn;
	double x0, x1, y0, y1;
	
	do{
		In = 0;
		Ipn = 0; 
		for(int i = 0;i < n;i++){
			x0 = a + i*(b-a)/n;
			x1 = x0 + (b-a)/n;
			y0 = f(x0);
			y1 = f(x1);
			In += y0 + y1;
		}
		In *= (b-a)/(2*n);
		for(int i = 0;i < q;i++){
			x0 = a + i*(b-a)/q;
			x1 = x0 + (b-a)/q;
			y0 = f(x0);
			y1 = f(x1);
			Ipn += y0 + y1;
		}
		Ipn *= (b-a)/(2*q);
		// cout << n << '\n';
		n *= p;
		q = n*p;
		// cout <<  (abs(In - Ipn)/(p*p + 1)) << '\n';
	}while((abs(In - Ipn)/(p*p + 1)) >= eps);
	return In; 
}

double hinhthang_case3(vector<double> &X, vector<double> &Y, double &eps){
	int q = X.size()-1;
	int p = 2;
	
	double In = 0, Ipn = 0;
	double eps;
	for(int i = 0;i < q;i++){
		Ipn += Y[i] + Y[i+1];
	}
	
	int k = 0;
	while(k < p){
		In += Y[k] + Y[k+p];
		k += p;
	}

	Ipn *= (X[1] - X[0])/2;
	In *= (X[p] - X[0])/2;

	eps = abs(In - Ipn)/(p*p +1); 
	return In; 
}
