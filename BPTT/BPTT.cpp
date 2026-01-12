// #include "../LinearEquationSystem/Gauss_Seiden.cpp"
#include "../Gauss-Jordan/Gauss-Jordan.cpp"

using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
vector<double> bptt_identify_function(vector<double> &X, vector<double> &Y);
vector<double> bptt_identify_function_case2(vector<double> &X, vector<double> &Y);
vector<double> bptt_identify_function_case3(vector<double> &X, vector<double> &Y);

double F1(double x);
double F2(double x);
double F3(double x);

//------------------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(){
	int n = 104;
    vector<double> X(n);
    vector<double> Y(n);
    for(int i = 0;i < n;i++){
    	cin >> X[i] >> Y[i];
    }
    
    vector<double> R = bptt_identify_function_case2(X,Y);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<double> bptt_identify_function(vector<double> &X, vector<double> &Y){
	int n = X.size() - 1;
	int m = 3;
	matrix P(n,m);

	for(int i = 0;i < n;i++){
		P.mat[i][0] = F1(X[i]);
	}
	for(int i = 0;i < n;i++){
		P.mat[i][1] = F2(X[i]);
	}
	for(int i = 0;i < n;i++){
		P.mat[i][2] = F3(X[i]);
	}
	cout << "P : \n";
	P.out();

	matrix M = P.t()*P;
	matrix Y_m(n,1);
	for(int i = 0;i < n;i++){
		Y_m.mat[i][0] = Y[i];
	}
	matrix Z = P.t()*Y_m;
	cout << "M : \n";
	M.out();
	cout << "Z : \n";
	Z.out();
	
	matrix A = gauss_jordan(M, Z);
	cout << "A : \n";
	A.out();
	vector<double> w = A.t().mat[0];
	// A.out();
	double p = 0; 
	double s,r;
	for(int i = 0;i < n;i++){
		s = w[0]*F1(X[i]) + w[1]*F2(X[i]) - Y[i];

		r = pow(s,2);
		p += r;
	}
	p = sqrt(p/n);
	cout << "sai so : " << p << '\n';

	return w;
}

vector<double> bptt_identify_function_case2(vector<double> &X, vector<double> &Y){
	int n = Y.size()-1;
	vector<int> K(0);
    vector<double> Y1(n);
    for(int i = 0;i < n;i++){
    	Y1[i] = log(abs(Y[i]));
    	if(Y[i] < 0){
    		K.push_back(i);
    	}
    	if(Y[i] == 0){
    		cout << "Zero in : " << i << '\n';
    		return X;
    	}
    }

    // for(auto i : Y1){
    // 	cout << i << ' ';
    // }

    vector<double> C = bptt_identify_function(X,Y1);
    C.back() = exp(C.back());

    for(auto i : C){
    	cout << i << ' ';
    }

   	return C;  
}

vector<double> bptt_identify_function_case3(vector<double> &X, vector<double> &Y){
	int n = Y.size()-1;
	vector<int> K_X(0);
	vector<int> K_Y(0);
	vector<double> X1(n);
    vector<double> Y1(n);
    for(int i = 0;i < n;i++){
    	X1[i] = abs(X[i]);
    	Y1[i] = log(abs(Y[i]));
    	if(Y[i] < 0){
    		K_Y.push_back(i);
    	}
    	if(Y[i] < 0){
    		K_X.push_back(i);
    	}
    	if(Y[i] == 0 || X[i] == 0){
    		cout << "Zero in : " << i << '\n';
    		return X;
    	}
    }

    // for(auto i : Y1){
    // 	cout << i << ' ';
    // }

    vector<double> C = bptt_identify_function(X1,Y1);
    C.back() = exp(C.back());

    for(auto i : C){
    	cout << i << ' ';
    }

   	return C;  
}

double F1(double x){
	return x*x;
}
double F2(double x){
	return 1/x;
}
double F3(double x){
	return 1;
}

