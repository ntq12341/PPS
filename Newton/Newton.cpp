#include "../Lagrange/Lagrange.cpp"
#include<cmath>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------
// thuật toán trích xuất khoảng nội suy tối ưu (cách đều)
vector<int> newton_extract_points(vector<double> &X, double x, int k);

// thuật toán tính tỉ sai phân
vector<double> newton_divided_diff(vector<double> &X,vector<double> &Y);

// thuật toán tính sai phân tiến
vector<double> newton_forward_diff(vector<double> &Y);

// thuật toán tính sai phân lùi
vector<double> newton_back_diff(vector<double> &Y);

// đa thức nội suy Newton tiến mốc bất kỳ
vector<double> newton_forward_non_equidistant(vector<double> &X,vector<double> &Y);

// đa thức nội suy Newton lùi mốc bất kỳ
vector<double> newton_back_non_equidistant(vector<double> &X,vector<double> &Y);

// đa thức nội suy Newton tiến mốc cách đều
vector<double> newton_forward_equidistant(vector<double> &Y);

// đa thức nội suy Newton lùi mốc cách đều
vector<double> newton_back_equidistant(vector<double> &Y);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// int main(){
//     vector<double> X = {0,1,2,3};
//     vector<double> Y = {1,4,9,16};
//     vector<double> S = newton_back_diff(Y);
//     for(auto i : S){
//         cout << i << ' ';
//     }
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<int> newton_extract_points(vector<double> &X, double x, int k){
    int n = X.size()-1; 
    double d = (x - X[0])/(X[1]-X[0]);
    // cout << d << '\n';
    vector<int> ans(0);
    int u,i1,i2,sp,ep;
    u = floor(d);
    if(d - round(d) < 1e-6){
        u = round(d);
    }
    // cout << u << '\n';
    if(k % 2 != 0){
        i1 = u - (k-1)/2;
        i2 = u + (k-1)/2;
    }else{
        i1 = u + 1 - k/2;
        i2 = u + k/2;
    }
    if(i1 < 0){
        sp = 0; 
        ep = k-1;
    }else if(i2 > n){
        sp = n-k+1;
        ep = n;
    }else{
        sp = i1; 
        ep = i2;
    }
    for(int j = sp;j <= ep;j++){
        ans.push_back(j);
    }
    return ans;
}

vector<double> newton_divided_diff(vector<double> &X,vector<double> &Y){
    int n = X.size()-1;
    int i = 1;
    vector<double> A = Y;
    vector<double> B;
    vector<double> ans;
    vector<vector<double>> Table; 
    while(i <= n){
        B.resize(n-i+1);
        for(int j = 0;j <= n-i;j++){
            B[j] = (A[j+1] - A[j])/(X[j+i] - X[j]);
        }
        ans.push_back(B[0]);
        Table.push_back(B);
        A = B;
        B.resize(0);
        i++;
    }
    // cout << "bang ty sai phan :" << '\n';
    // for(auto v : Table){
    //     for(auto i : v){
    //         cout << i << ' ';
    //     }
    //     cout << '\n';
    // }
    return ans; 
}

vector<double> newton_forward_diff(vector<double> &Y){
    int n = Y.size()-1;
    int i = 1;
    vector<double> A = Y;
    vector<double> B;
    vector<double> ans;
    vector<vector<double>> Table; 
    while(i <= n){
        B.resize(n-i+1);
        for(int j = 0;j <= n-i;j++){
            B[j] = (A[j+1] - A[j]);
        }
        ans.push_back(B[0]);
        Table.push_back(B);
        A = B;
        B.resize(0);
        i++;
    }
    // cout << "bang sai phan tien :" << '\n';
    // for(auto v : Table){
    //     for(auto i : v){
    //         cout << i << ' ';
    //     }
    //     cout << '\n';
    // }
    return ans;
}

vector<double> newton_back_diff(vector<double> &Y){
    int n = Y.size()-1;
    int i = 1;
    vector<double> A = Y;
    vector<double> B;
    vector<double> ans;
    vector<vector<double>> Table; 
    while(i <= n){
        B.resize(n-i+1);
        for(int j = n;j >= i;j--){
            B[j-i] = (A[j-i+1] - A[j-i]);
        }
        ans.push_back(B.back());
        Table.push_back(B);
        A = B;
        B.resize(0);
        i++;
    }
    // cout << "bang sai phan lui :" << '\n';
    // for(auto v : Table){
    //     for(auto i : v){
    //         cout << i << ' ';
    //     }
    //     cout << '\n';
    // }
    return ans;
}

vector<double> newton_forward_non_equidistant(vector<double> &X,vector<double> &Y){
    int n = X.size()-1;
    vector<double> F = newton_divided_diff(X,Y);
    int i = 1;
    vector<double> P(n+1);
    P[0] = Y[0];
    vector<double> w;
    vector<double> Omega;
    while(i <= n){
        // cout << "i : " << i << '\n';
        // cout << "P : \n";
        // for(auto z : P){
        //     cout << z << ' '; 
        // } 
        // cout << '\n';
        w.resize(0);
        for(int j = 0;j < i;j++){
            w.push_back(X[i]);
        }
        Omega = horner_omega(w);
        for(int j = 0;j <= i;j++){
            P[j] += F[i-1]*Omega[j];
        }
        i++;
    }
    // cout << "P : \n";
    // for(auto z : P){
    //     cout << z << ' '; 
    // } 
    // cout << '\n';
    return P;
}

vector<double> newton_back_non_equidistant(vector<double> &X,vector<double> &Y){
    int n = X.size()-1;
    vector<double> X_i(0);
    vector<double> Y_i(0);
    for(int j = n;j >= 0;j--){
        X_i.push_back(X[j]);
        Y_i.push_back(Y[j]);
    }
    vector<double> F = newton_divided_diff(X_i,Y_i);
    int i = 1;
    vector<double> P(n+1);
    P[0] = Y_i[0];
    vector<double> w;
    vector<double> Omega;
    while(i <= n){
        cout << "i : " << i << '\n';
        cout << "P : \n";
        for(auto z : P){
            cout << z << ' '; 
        } 
        cout << '\n';
        w.resize(0);
        for(int j = 0;j < i;j++){
            w.push_back(X_i[i]);
        }
        Omega = horner_omega(w);
        for(int j = 0;j <= i;j++){
            P[j] += F[i-1]*Omega[j];
        }
        i++;
    }
    cout << "i : " << i << '\n';
    cout << "P : \n";
    for(auto z : P){
        cout << z << ' '; 
    } 
    cout << '\n';
    return P;
}

vector<double> newton_forward_equidistant(vector<double> &Y){
    int n = Y.size()-1;
    vector<double> D = newton_forward_diff(Y);
    int i = 1;
    vector<double> P(n+1);
    P[0] = Y[0];
    vector<double> w;
    vector<double> Omega;
    while(i <= n){
        w.resize(0);
        for(int j = 0;j < i;j++){
            w.push_back(j);
        }
        Omega = horner_omega(w);
        for(int j = 0;j <= i;j++){
            P[j] += (D[i-1]*Omega[j])/factorial(i);
        }
        i++;
    }
    return P; 
}

vector<double> newton_back_equidistant(vector<double> &Y){
    int n = Y.size()-1;
    vector<double> D = newton_back_diff(Y);
    int i = 1;
    vector<double> P(n+1);
    P[0] = Y[n];
    vector<double> w;
    vector<double> Omega;
    while(i <= n){
        w.resize(0);
        for(int j = 0;j < i;j++){
            w.push_back(-j);
        }
        Omega = horner_omega(w);
        for(int j = 0;j <= i;j++){
            P[j] += (D[i-1]*Omega[j])/factorial(i);
        }
        i++;
    }
    return P; 
}