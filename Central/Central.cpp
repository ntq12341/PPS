#include "../Newton/Newton.cpp"
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------

// nội suy trung tâm Gauss I 
vector<double> central_gaussI(vector<double> &Y);
// nội suy trung tâm Gauss II 
vector<double> central_gaussII(vector<double> &Y);
// nội suy trung tâm Sterling
vector<double> central_sterling(vector<double> &Y);
// nội suy trung tâm Bessel 
vector<double> central_bessel(vector<double> &Y);
// bảng sai phân tiến 
vector<vector<double>> central_forward_diff(vector<double> &Y);

//------------------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// int main(){
//     vector<double> Y = {4,1,0,1,4,9};
//     // vector<vector<double>> S = central_forward_diff(Y);
//     // for(auto a : S){
//     //     for(auto i : a){
//     //         cout << i << ' ';
//     //     }
//     //     cout << '\n';
//     // }
//     vector<double> S = central_bessel(Y);
//     for(auto i : S){
//         cout << i << ' ';
//     }
//     cout << '\n';
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<double> central_gaussI(vector<double> &Y){
    int n = (Y.size()-1)/2;
    vector<vector<double>> Delta = central_forward_diff(Y);
    // for (auto i : Delta){
    //     for(auto j : i){
    //         cout << j << ' ';
    //     }
    //     cout << '\n';
    // }
    vector<double> P(2*n+1, 0);
    P[0] = Y[n];
    P[1] = Delta[0][n];
    int i = 1;
    double a_0, a_1, u;
    vector<double> T;
    while(i <= n){

        // cout << "i : " << i << '\n';

        a_0 = Delta[2*i-1][n-i]/factorial(2*i);

        // cout << "a_2i : " << a_0 << '\n';

        T = vector<double>{1};
        for(int j = 0; j <= i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 1;j < i;j++){
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_0*T[j]; 
        }

        if(i == n){
            break;
        }

        a_1 = Delta[2*i][n-i]/factorial(2*i+1);

        // cout << "a_2i+1 : " << a_1 << '\n';

        T = vector<double>{1};
        for(int j = 0; j <= i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 1;j <= i;j++){
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_1*T[j];
        }

        // cout << "P : \n" ;
        // for (auto k : P){
        //     cout << k << ' ';
        // }
        // cout << '\n';


        i = i+1;
    }

    // cout << "P : \n"; 
    // for (auto k : P){
    //     cout << k << ' ';
    // }
    // cout << '\n';

    return P; 
}

vector<double> central_gaussII(vector<double> &Y){
    int n = (Y.size()-1)/2;
    vector<vector<double>> Delta = central_forward_diff(Y);
    vector<double> P(2*n+1, 0);
    P[0] = Y[n];
    int i = 1;
    double a_0, a_1, u;
    vector<double> T;
    while(i <= n){
        a_0 = Delta[2*i-1][n-i]/factorial(2*i);
        T = vector<double>{1};
        for(int j = 0; j <= i;j++){
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 1;j < i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_0*T[j]; 
        }

        a_1 = Delta[2*i-2][n-i]/factorial(2*i-1);
        T = vector<double>{1};
        for(int j = 0; j < i;j++){
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 1;j < i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_1*T[j];
        }

        i = i+1;
    }
    return P; 
}

vector<double> central_sterling(vector<double> &Y){
    int n = (Y.size()-1)/2;
    vector<vector<double>> Delta = central_forward_diff(Y);
    vector<double> P(2*n+1, 0);
    P[0] = Y[n];
    int i = 1;
    double a_0, a_1, u;
    vector<double> T;
    while(i <= n){ 
        a_0 = Delta[2*i-1][n-i]/(factorial(2*i));

        T = vector<double>{0,0,1};
        for(int j = 1;j < i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_0*T[j];
        }

        a_1 = (Delta[2*i-2][n-i] + Delta[2*i-2][n-i+1])/(2*factorial(2*i-1));

        T = vector<double>{0,1};
        for(int j = 1;j < i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
            u = -j;
            T = horner_polynominal_multiply(T,u);
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_1*T[j];
        }

        i = i+1;
    }
    return P; 
}

vector<double> central_bessel(vector<double> &Y){
    int n = (Y.size()-2)/2;
    vector<vector<double>> Delta = central_forward_diff(Y);
    vector<double> P(2*n+2, 0);
    P[0] = (Y[n]+ Y[n+1]-Delta[0][n])/2;
    P[1] = Delta[0][n];
    int i = 1;
    double a_0, a_1, u;
    vector<double> T;
    while(i <= n){ 

        cout << "i : " << i << '\n';
        
        a_0 = (Delta[2*i-1][n-i] + Delta[2*i-1][n-i+1])/(2*factorial(2*i));

        cout << "a_2i : " << a_0 << '\n';

        T = vector<double>{1};
        for(int j = 0;j <= i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
            if(j >= 1 && j < i){
                u = -j;
                T = horner_polynominal_multiply(T,u);
            }
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_0*T[j];
        }

        a_1 = Delta[2*i][n-i]/(factorial(2*i+1));

        cout << "a_2i+1 : " << a_1 << '\n';

        T = vector<double>{-0.5,1};
        for(int j = 0;j <= i;j++){
            u = j;
            T = horner_polynominal_multiply(T,u);
            if(j >= 1 && j < i){
                u = -j;
                T = horner_polynominal_multiply(T,u);
            }
        }
        for(int j = 0; j < T.size(); j++){
            P[j] += a_1*T[j];
        }

        cout << "P : \n"; 
        for (auto k : P){
            cout << k << ' ';
        }
        cout << '\n';

        i = i+1;
    }
    cout << "P : \n"; 
    for (auto k : P){
        cout << k << ' ';
    }
    cout << '\n';
    return P; 
}

vector<vector<double>> central_forward_diff(vector<double> &Y){
    int n = Y.size()-1;
    int i = 1;
    vector<double> A = Y;
    vector<double> B;
    vector<vector<double>> ans;
    while(i <= n){
        B.resize(n-i+1);
        for(int j = 0;j <= n-i;j++){
            B[j] = (A[j+1] - A[j]);
        }
        ans.push_back(B);
        A = B;
        B.resize(0);
        i++;
    }
    return ans;
}