#include "../Horner/Horner.cpp"
#include <math.h>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
// xấp xỉ giá trị đạo hàm 
double dao_ham_calculate(vector<double> &X,vector<double> &Y,double x){
	int p = Y.size()-1; 
	double h = X[1]-X[0];
	double ov_t = (x-X[0])/h;
	vector<double> A(p+1,0);
	vector<double> Tmp(0);
	double ans = 0;
	double c = 0;
	for(int i = 0;i <= p;i++){
		Tmp.push_back(i);  
	}
	vector<double> Omega = horner_omega(Tmp); 
	for(int i = 0;i <= p;i++){
		vector<double> T = horner_polynominal_divide(Omega,i,c);
		vector<double> T1(p);  
		for(int j = 0;j < p;j++){
			T1[j] += T[j+1]*(j+1);
		}
		// for(auto i : T1){
		// 	cout << i << ' ';
		// }
		// cout << '\n'; 
		c = pow(-1,p-i)/(factorial(i)*factorial(p-i));
		cout << c << '\n'; 
		A[i] = c*horner_polynominal_value(T1,ov_t); 

		ans += A[i]*Y[i];
		c = 0;
	}
	ans /= h; 
	return ans; 
}

