#include <iostream>
#include <vector>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------

// tính giá trị đa thức P tại c
double horner_polynominal_value(vector<double> &P, double c);

// nhân đa thức P với x-c
vector<double> horner_polynominal_multiply(vector<double> &P, double c);

// chia đa thức P cho x-c lấy số dư r
vector<double> horner_polynominal_divide(vector<double> &P, double c, double &r);

// tính đạo hàm cấp k của đa thức P tại c 
double horner_polynominal_derivation(vector<double> &P, int k, double c);

// tính các hệ số của đa thức omega của tập P gồm n điểm
vector<double> horner_omega(vector<double> &P);
//------------------------------------------------------------------------

// tính giai thừa
ll factorial(int n);
// thêm 1 phần tử vào đầu mảng 
void push_head(vector<double> &v, double & k);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// int main(){
//     vector<double> P = {1,2,3};
//     vector<double> S = horner_omega(P);
//     for(auto i : S){
//         cout << i << ' ';
//     }
//    	return 0;
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double horner_polynominal_value(vector<double> &P, double c){
    int n = P.size()-1;
    int k = n;
    double b = P[n];
    while(k > 0){
        // cout << "k : " << k << '\n';
        // cout << "b : " << b << '\n';
        k--;
        b = P[k] + b*c;
        // cout << "b : " << b << '\n';
    }
    // cout << "k : " << k << '\n';
    // cout << "b : " << b << '\n';
    return b; 
}

vector<double> horner_polynominal_multiply(vector<double> &P, double c){
    int n = P.size() - 1;
    int k = 0;
    vector<double> ans;
    ans.push_back(-c*P[k++]);
    while(k <= n){
        ans.push_back(P[k-1] - c*P[k++]);
    }
    ans.push_back(P[k-1]);
    return ans; 
}

vector<double> horner_polynominal_divide(vector<double> &P, double c, double &r){
    r = 0;
    int n = P.size()-1;
    int k = n;
    vector<double> res;
    res.push_back(P[k]);
    while(k > 1){
        // cout << k << " -res : \n" ;
        // for(auto i : res){
        //     cout << i << ' ';
        // }
        // cout << '\n';
        k--;
        r = c*res[0] + P[k];
        push_head(res, r);
    }
    r = c*res[0] + P[0];
    // cout << "thuong : \n";
    // for(auto i : res){
    //     cout << i << ' ';
    // }
    // cout << '\n';
    // cout << "du : " << r << '\n';
    return res;
}

double horner_polynominal_derivation(vector<double> &P, int k, double c){
    int n = P.size()-1;
    int j = 0;
    double r = 0; 
    vector<double> Q; 
    vector<double> D = P;
    vector<double> R(0); 
    while(j <= n){
        // cout << "J : " << j << "\n D : \n" ;
        // for(auto i : D){
        //     cout << i << ' ';
        // }
        // cout << "\n Q : \n";
        Q = horner_polynominal_divide(D,c,r);
        // for(auto i : Q){
        //     cout << i << ' ';
        // }
        D = Q; 
        // cout << "\n R : \n";
        R.push_back(r);
        // for(auto i : R){
        //     cout << i << ' ';
        // }
        // cout << '\n';
        j++; 
    }
    return R[k]*factorial(k);

}

vector<double> horner_omega(vector<double> &P){
    int n = P.size()-1;
    vector<double> ans = {-P[0], 1};
    vector<double> R;
    int k = 1;
    while(k <= n){
        // cout << "k : " << k << " - x : " << P[k-1] << "\n A :";
        // for (auto i : ans){
        //     cout << i << ' ';
        // }
        // cout << '\n';
        R = horner_polynominal_multiply(ans,P[k++]); 
        ans = R;
    }
    // cout << "k : " << k << " - x : " << P[k-1] << "\n A :";
    // for (auto i : ans){
    //     cout << i << ' ';
    // }
    // cout << '\n';
    return ans;
}
//------------------------------------------------------------------------
void push_head(vector<double> &v, double & k){
    v.push_back(0);
    for(int i = v.size()-2;i >= 0;i--){
        v[i+1] = v[i];
    }
    v[0] = k;
}

ll factorial(int n) {
    if (n < 0) return -1; // Không tồn tại giai thừa cho số âm
    if (n == 0 || n == 1) return 1;
    
    ll result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}