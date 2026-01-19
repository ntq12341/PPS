#include "../Horner/Horner.cpp"
#include "../matrix.h" 
#include <math.h>
using namespace std;

double f(double x,double y){
	return pow(y,3)*x*sin(x+y); 
}

vector<vector<double>> RK1(double y0, double x0, double xn, double h){
	double r1 = 1;
	double k1_n = 0;

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	double y_tmp = y0;
	double x_tmp = x0; 

	while(x_tmp < xn){
		k1_n = h*f(x_tmp, y_tmp);
		y_tmp += r1*k1_n; 
		x_tmp += h;

		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> RK2(double y0, double x0, double xn, double h){
	double r1 = 0.5, r2 = 0.5;
	double k1_n, k2_n;
	double a2 = 1,b11 = 1;

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	double y_tmp = y0;
	double x_tmp = x0; 

	while(x_tmp < xn){
		k1_n = h*f(x_tmp, y_tmp);
		k2_n = h*f(x_tmp + a2*h, y_tmp + b11*k1_n);
		y_tmp += r1*k1_n + r2*k2_n; 
		x_tmp += h;

		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> RK3(double y0, double x0, double xn, double h){
	double r1 = 1.0/6, r2 = 2.0/3, r3 = 1.0/6;
	double k1_n, k2_n, k3_n;
	double a2 = 0.5,a3 = 1;
	double b11 = 0.5, b21 = -1, b22 = 2;

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	double y_tmp = y0;
	double x_tmp = x0; 

	while(x_tmp < xn){ 
		k1_n = h*f(x_tmp, y_tmp);
		k2_n = h*f(x_tmp + a2*h, y_tmp + b11*k1_n);
		k3_n = h*f(x_tmp + a3*h, y_tmp + b21*k1_n + b22*k2_n);
		y_tmp += r1*k1_n + r2*k2_n + r3*k3_n; 
		x_tmp += h;
		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> RK4(double y0, double x0, double xn, double h){
	double r1 = 1.0/6, r2 = 1.0/3, r3 = 1.0/3, r4 = 1.0/6;
	double k1_n, k2_n, k3_n, k4_n;
	double a2 = 0.5, a3 = 0.5, a4 = 1; 
	double b11 = 0.5, b21 = 0, b22 = 0.5, b31 = 0, b32 = 0, b33 = 1;

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	double y_tmp = y0;
	double x_tmp = x0; 

	while(x_tmp < xn){
		k1_n = h*f(x_tmp, y_tmp);
		k2_n = h*f(x_tmp + a2*h, y_tmp + b11*k1_n);
		k3_n = h*f(x_tmp + a3*h, y_tmp + b21*k1_n + b22*k2_n);
		k4_n = h*f(x_tmp + a4*h, y_tmp + b31*k1_n + b32*k2_n + b33*k3_n);
		y_tmp += r1*k1_n + r2*k2_n + r3*k3_n + r4*k4_n; 
		x_tmp += h;

		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}