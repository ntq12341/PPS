#include "../RK/RK.cpp"
using namespace std;

// double f(double x,double y){
// 	return pow(y,3)*x*sin(x+y); 
// }

vector<vector<double>> Adam2(double y0, double x0, double xn, double h){
	vector<vector<double>> RK = RK4(y0, x0, xn, h);

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	vector<double> F(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	F.push_back(f(x0,y0));

	double y_tmp = y0, y_bef;
	double x_tmp = x0; 

	int i;

	while(x_tmp < xn){
		x_tmp += h;
		i = (x_tmp - x0)/(double)h; 
		if(i < 2){
			F.push_back(f(x_tmp, RK[1][i]));
			continue;
		} 
		y_bef = y_tmp;
		// AB
		y_tmp = y_bef + (h/2.0)*(3*F[i-1] - F[i-2]);
		// AM
		for(int k = 0;k < 10;k++){
			y_tmp = y_bef + (h/12.0)*(5*f(x_tmp,y_tmp) - 8*F[i-1] - F[i-2]);
		}

		F.push_back(f(x_tmp,y_tmp));
		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> Adam3(double y0, double x0, double xn, double h){
	vector<vector<double>> RK = RK4(y0, x0, xn, h);

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	vector<double> F(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	F.push_back(f(x0,y0));

	double y_tmp = y0, y_bef;
	double x_tmp = x0; 

	int i;

	while(x_tmp < xn){
		x_tmp += h;
		i = (x_tmp - x0)/(double)h; 
		if(i < 3){
			F.push_back(f(x_tmp, RK[1][i]));
			continue;
		} 
		y_bef = y_tmp;
		// AB
		y_tmp = y_bef + (h/12.0)*(23*F[i-1] - 16*F[i-2] + 5*F[i-3]);
		// AM
		for(int k = 0;k < 10;k++){
			y_tmp = y_bef + (h/24.0)*(9*f(x_tmp,y_tmp) + 19*F[i-1] - 5*F[i-2] + F[i-3]);
		}

		F.push_back(f(x_tmp,y_tmp));
		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> Adam4(double y0, double x0, double xn, double h){
	vector<vector<double>> RK = RK4(y0, x0, xn, h);

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	vector<double> F(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	F.push_back(f(x0,y0));

	double y_tmp = y0, y_bef;
	double x_tmp = x0; 

	int i;

	while(x_tmp < xn){
		x_tmp += h;
		i = (x_tmp - x0)/(double)h; 
		if(i < 4){
			F.push_back(f(x_tmp, RK[1][i]));
			continue;
		} 
		y_bef = y_tmp;
		// AB
		y_tmp = y_bef + (h/24.0)*(55*F[i-1] - 59*F[i-2] + 37*F[i-3] - 9*F[i-4]);
		// AM
		for(int k = 0;k < 10;k++){
			y_tmp = y_bef + (h/720.0)*(251*f(x_tmp,y_tmp) + 646*F[i-1] - 264*F[i-2] + 106*F[i-3] - 19*F[i-4]);
		}

		F.push_back(f(x_tmp,y_tmp));
		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}

vector<vector<double>> Adam5(double y0, double x0, double xn, double h){
	vector<vector<double>> RK = RK4(y0, x0, xn, h);

	vector<double> X_ans(0);
	vector<double> Y_ans(0);

	vector<double> F(0);

	Y_ans.push_back(y0);
	X_ans.push_back(x0); 

	F.push_back(f(x0,y0));

	double y_tmp = y0, y_bef;
	double x_tmp = x0; 

	int i;

	while(x_tmp < xn){
		x_tmp += h;
		i = (x_tmp - x0)/(double)h; 
		if(i < 5){
			F.push_back(f(x_tmp, RK[1][i]));
			continue;
		} 
		y_bef = y_tmp;
		// AB
		y_tmp = y_bef + (h/720.0)*(1901*F[i-1] - 2774*F[i-2] + 2616*F[i-3] - 1274*F[i-4] + 251*F[i-5]);
		// AM
		for(int k = 0;k < 10;k++){
			y_tmp = y_bef + (h/1440.0)*(475*f(x_tmp,y_tmp) + 1427*F[i-1] - 798*F[i-2] + 482*F[i-3] - 173*F[i-4] + 27*F[i-5]);
		}

		F.push_back(f(x_tmp,y_tmp));
		X_ans.push_back(x_tmp);
		Y_ans.push_back(y_tmp);
	}
	vector<vector<double>> ans = {X_ans,Y_ans};
	return ans;
}
