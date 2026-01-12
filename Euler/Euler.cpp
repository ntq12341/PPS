#include "../Horner/Horner.cpp"
#include <math.h>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
double F(double x,double y){
	return pow(y,3)*x*sin(x+y); 
}

double dF_dy(double x, double y){
    return x * ( 3*pow(y,2)*sin(x+y) + pow(y,3)*cos(x+y) );
}

vector<vector<double>> Euler_Hien(double y0, double x0, double xn, double h){
	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	double y_tmp = y0;
	double x_tmp = x0; 

	while(x_tmp < xn){
		y_tmp += F(x_tmp, y_tmp)*h; 
		x_tmp += h;

		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans; 
}

vector<vector<double>> Euler_An(
    double y0, double x0, double xn, double h,
    int max_iter = 200, double tol = 1e-10
){
    vector<double> X_ans, Y_ans;

    double x = x0;
    double y = y0;

    X_ans.push_back(x);
    Y_ans.push_back(y);

    int N = (int)((xn - x0)/h);

    for(int n = 0; n <= N; n++){
        double x_next = x + h;

        // initial guess: Euler hiện
        double y_guess = y + h * F(x, y);

        // Newton iteration
        for(int k = 0; k < max_iter; k++){
            double g  = y_guess - y - h * F(x_next, y_guess);
            double dg = 1.0 - h * dF_dy(x_next, y_guess);

            double y_new = y_guess - g / dg;

            if (fabs(y_new - y_guess) < tol){
                y_guess = y_new;
                break;
            }
            y_guess = y_new;
        }

        y = y_guess;
        x = x_next;

        X_ans.push_back(x);
        Y_ans.push_back(y);
    }

    return {X_ans, Y_ans};
}

vector<vector<double>> Euler_CaiTien(
    double y0, double x0, double xn, double h,
    int max_iter = 200, double tol = 1e-10
){
    vector<double> X, Y;

    double x = x0;
    double y = y0;

    X.push_back(x);
    Y.push_back(y);

    int N = (int)((xn - x0) / h);

    for(int n = 0; n <= N; n++){
        double x_next = x + h;

        // predictor: Euler hiện
        double y_guess = y + h * F(x, y);

        // Newton–Raphson
        for(int k = 0; k < max_iter; k++){
            double g =
                y_guess - y
                - 0.5*h*( F(x, y) + F(x_next, y_guess) );

            double dg =
                1.0 - 0.5*h * dF_dy(x_next, y_guess);

            double y_new = y_guess - g / dg;

            if (fabs(y_new - y_guess) < tol){
                y_guess = y_new;
                break;
            }
            y_guess = y_new;
        }

        y = y_guess;
        x = x_next;

        X.push_back(x);
        Y.push_back(y);
    }

    return {X, Y};
}
