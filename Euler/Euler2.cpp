#include "../Horner/Horner.cpp"
#include "../matrix.h" 
#include <math.h>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
double F1(double x,double y, double z){
	return 0.4*x*(1-x/20) + 0.4*y - 0.3*x*z; 
}

double F2(double x,double y, double z){
	return 0.7*y*(1-y/25) - 0.4*y - 0.4*y*z; 
}

double F3(double x, double y, double z){
	return -0.3*z + 0.35*(x+y)*z;
}

double dF1_x(double x,double y, double z){
	return 0.4 - 0.04*x - 0.3*z; 
}

double dF1_y(double x,double y, double z){
	return 0.4;
}

double dF1_z(double x,double y, double z){
	return -0.3*x;
}

double dF2_x(double x,double y, double z){
	return 0;
}

double dF2_y(double x,double y, double z){
	return 0.3 - (1.4*y)/25 - 0.4*z;
}

double dF2_z(double x,double y, double z){
	return -0.4*y;
}

double dF3_x(double x,double y, double z){
	return 0.35*z;
}

double dF3_y(double x,double y, double z){
	return 0.35*z;
}

double dF3_z(double x,double y, double z){
	return -0.3 + 0.35*(x+y);
}

matrix J(double x, double y, double z, double h){
	matrix ans(3);

	ans.mat[0][0] = 1 - h*dF1_x(x,y,z);
	ans.mat[0][1] = -h*dF1_y(x,y,z);
	ans.mat[0][2] = -h*dF1_z(x,y,z);
	ans.mat[1][0] = -h*dF2_x(x,y,z);
	ans.mat[1][1] = 1 - h*dF2_y(x,y,z);
	ans.mat[1][2] = -h*dF2_z(x,y,z);
	ans.mat[2][0] = -h*dF3_x(x,y,z);
	ans.mat[2][1] = -h*dF3_y(x,y,z);
	ans.mat[2][2] = 1 - h*dF3_z(x,y,z);

	return ans; 
}

vector<vector<double>> Euler_Hien_Hept(
    double x0, double y0, double z0, double t0, double tn, double h){
    vector<double> T_ans, X_ans, Y_ans, Z_ans;

    double t = t0;
    double x = x0;
    double y = y0;
    double z = z0;

    T_ans.push_back(t);
    X_ans.push_back(x);
    Y_ans.push_back(y);
    Z_ans.push_back(z);

    int N = (int)((tn - t0)/h);

    for(int n = 0; n <= N; n++){
        double t_next = t + h;

        // initial guess: Euler hiện
        double x_tmp = x + h * F1(x, y,z);
        double y_tmp = y + h * F2(x, y,z);
        double z_tmp = z + h * F3(x, y,z);

        x = x_tmp;
        y = y_tmp;
        z = z_tmp;

        X_ans.push_back(x);
        Y_ans.push_back(y);
        Z_ans.push_back(z);
        T_ans.push_back(t);
    }

    return {T_ans,X_ans,Y_ans,Z_ans};
}

vector<vector<double>> Euler_An_Hept(
    double x0, double y0, double z0, double t0, double tn, double h,
    int max_iter = 3, double tol = 1e-10
){
    vector<double> T_ans, X_ans, Y_ans, Z_ans;

    double t = t0;
    double x = x0;
    double y = y0;
    double z = z0;

    T_ans.push_back(t);
    X_ans.push_back(x);
    Y_ans.push_back(y);
    Z_ans.push_back(z);

    int N = (int)((tn - t0)/h);

    for(int n = 0; n <= N; n++){
        double t_next = t + h;

        // initial guess: Euler hiện
        double x_guess = x + h * F1(x, y,z);
        double y_guess = y + h * F2(x, y,z);
        double z_guess = z + h * F3(x, y,z);
        matrix guess(3,1);
        guess.mat[0][0] = x_guess;
        guess.mat[1][0] = y_guess;
        guess.mat[2][0] = z_guess;

        matrix J_i = J(x_guess,y_guess,z_guess,h).inverse(); 

        // Newton iteration
        for(int k = 0; k < max_iter; k++){
            double g_x  = x_guess - x - h * F1(x_guess, y_guess, z_guess);
            double g_y  = y_guess - y - h * F2(x_guess, y_guess, z_guess);
            double g_z  = z_guess - z - h * F3(x_guess, y_guess, z_guess);
            matrix g(3,1);
            g.mat[0][0] = g_x;
            g.mat[1][0] = g_y;
            g.mat[2][0] = g_z;

            matrix next = guess - J_i*g;
            guess = next; 
        }

        x = x_guess;
        y = y_guess;
        z = z_guess;
        t = t_next;

        X_ans.push_back(x);
        Y_ans.push_back(y);
        Z_ans.push_back(z);
        T_ans.push_back(t);
    }

    return {T_ans,X_ans,Y_ans,Z_ans};
}

vector<vector<double>> Euler_CaiTien_Hept(
    double x0, double y0, double z0, double t0, double tn, double h,
    int max_iter = 3, double tol = 1e-10
){
    vector<double> T_ans, X_ans, Y_ans, Z_ans;

    double t = t0;
    double x = x0;
    double y = y0;
    double z = z0;

    T_ans.push_back(t);
    X_ans.push_back(x);
    Y_ans.push_back(y);
    Z_ans.push_back(z);

    int N = (int)((tn - t0)/h);

    for(int n = 0; n <= N; n++){
        double t_next = t + h;

        // initial guess: Euler hiện
        double x_guess = x + h * F1(x, y,z);
        double y_guess = y + h * F2(x, y,z);
        double z_guess = z + h * F3(x, y,z);
        matrix guess(3,1);
        guess.mat[0][0] = x_guess;
        guess.mat[1][0] = y_guess;
        guess.mat[2][0] = z_guess;

        matrix J_i = J(x_guess,y_guess,z_guess,h).inverse(); 

        // Newton iteration
        for(int k = 0; k < max_iter; k++){
            double g_x  = x_guess - x - (h/2) * (F1(x_guess, y_guess, z_guess) + F1(x,y,z));
            double g_y  = y_guess - y - (h/2) * (F2(x_guess, y_guess, z_guess) + F2(x,y,z));
            double g_z  = z_guess - z - (h/2) * (F3(x_guess, y_guess, z_guess) + F3(x,y,z));
            matrix g(3,1);
            g.mat[0][0] = g_x;
            g.mat[1][0] = g_y;
            g.mat[2][0] = g_z;

            matrix next = guess - J_i*g;
            guess = next; 
        }

        x = x_guess;
        y = y_guess;
        z = z_guess;
        t = t_next;

        X_ans.push_back(x);
        Y_ans.push_back(y);
        Z_ans.push_back(z);
        T_ans.push_back(t);
    }

    return {T_ans,X_ans,Y_ans,Z_ans};
}