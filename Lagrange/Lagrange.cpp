// https://www.dcode.fr/lagrange-interpolating-polynomial
//________________________________________________________
#include "../Horner/Horner.cpp"
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
// nội suy Lagrange mốc bất kỳ
vector<double> lagrange_non_equidistant(vector<double> &X,vector<double> &Y);

// nội suy Lagrange mốc cách đều
vector<double> lagrange_equidistant(vector<double> &Y);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// int main(){
//     vector<double> X = {1,2,3,4};
//     vector<double> Y = {4,5,6,7};
//     vector<double> S = lagrange_equidistant(Y);
//     for(auto i : S){
//     	cout << i << ' ';
//     }
//     cout << '\n';
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<double> lagrange_non_equidistant(vector<double> &X,vector<double> &Y){
	int n = X.size()-1;
	// for (auto i : X){
    //     cout << i << ' ' ;
    // }
    // cout << '\n';
	// cout << "n : " << n << '\n' ;
	vector<double> Omega = horner_omega(X);
	// cout << "omega : \n";
	// for (auto i : Omega){
    //     cout << i << ' ' ;
    // }
    // cout << '\n';
	int i = 0;
	vector<double> ans(n+1,0);
	vector<double> P;
	double tmp; 
	while(i <= n){
		// cout << "i : " << i << '\n';
		// cout << X[i] << '\n';
		P = horner_polynominal_divide(Omega,X[i],tmp);
		// cout << "P : \n";
		// for (auto i : P){
	    //     cout << i << ' ';
	    // }
	    // cout << '\n';
		tmp = Y[i];
		for(int j = 0;j <= n;j++){
			if(j != i){
				tmp /= (X[i] - X[j]);
			}
		}
		// cout << "R : " << tmp << '\n';
		for(int j = 0;j <= n;j++){
			ans[j] += P[j]*tmp;
		}
		// cout << "A : \n";
		// for (auto i : ans){
	    //     cout << i << ' ';
	    // }
	    // cout << '\n';
		i++;
	}
	// cout << "A : \n";
	// for (auto i : ans){
    //     cout << i << ' ';
    // }
    // cout << '\n';
	return ans; 
}

vector<double> lagrange_equidistant(vector<double> &Y){
	int n = Y.size()-1;
	vector<double> P;
	for(int j = 0;j <= n;j++){
		P.push_back(j);
	}
	double i = 0;
	vector<double> ans(n+1,0);
	vector<double> Omega = horner_omega(P);
	double tmp; 
	while(i <= n){
		P = horner_polynominal_divide(Omega,i,tmp);
		tmp = Y[i];
		for(int j = 0;j <= n;j++){
			if(j != i){
				tmp /= (i - j);
			}
		}
		for(int j = 0;j <= n;j++){
			ans[j] += P[j]*tmp;
		}
		i++;
	}
	return ans; 
}