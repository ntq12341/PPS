#include "../Central/Central.cpp"
#include <cmath>
using namespace std;

typedef long long ll;
//------------------------------------------------------------------------

// khoảng cách ly nghiệm 
vector<double> inverse_solution_range(vector<double> &Y,double y);
// xác định các khoảng đơn điệu 
vector<double> inverse_monotonic_range(vector<double> &Y);
// hàm ngược
vector<double> inverse_function(vector<double> &X,vector<double> &Y,double y);
// lặp newton tiến 
vector<double> inverse_newton_forward(vector<double> &X,vector<double> &Y,double y,double error);
// lặp newton lùi
vector<double> inverse_newton_back(vector<double> &X,vector<double> &Y,double y,double error);
// bảng sai phân tiến 
vector<vector<double>> inverse_forward_diff(vector<double> &Y);
// bảng sai phân lùi
vector<vector<double>> inverse_back_diff(vector<double> &Y);

//------------------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// int main(){
//     vector<double> X(0);
//     vector<double> Y(0);
//     double d;
//     for(int i = 0;i < 141;i++){
//         cin >> d;
//         X.push_back(d);
//     }
//     for(int i = 0;i < 141;i++){
//         cin >> d;
//         Y.push_back(d);
//     }
//     // for(auto i : X){
//     //     cout << i << ' ';
//     // }
//     // cout << '\n';
//     // for(auto i : Y){
//     //     cout << i << ' ';
//     // }
//     // cout << '\n';
//     vector<double> S = inverse_newton_forward(X,Y,3.15,1e-3);
//     for(auto i : S){
//         cout << i << ' ';
//     }
// }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<double> inverse_solution_range(vector<double> &Y,double y){
    vector<double> R(0);
    for(int i = 0;i < Y.size()-1;i++){
        if(Y[i] <= y && y <= Y[i+1]){
            R.push_back(i);
        }else if(Y[i] >= y && y >= Y[i+1]){
            R.push_back(i);
        }
    }
    return R;
}

vector<double> inverse_monotonic_range(vector<double> &Y){
    vector<double> R = {0};
    vector<double> Z(0);
    for(int i = 0;i < Y.size()-1;i++){
        Z.push_back( (Y[i+1]-Y[i]) / abs(Y[i]-Y[i+1]) );
    }
    for(int i = 1;i < Z.size();i++){
        if(Z[i]*Z[i-1] == -1){
            R.push_back(i);
        }
    }
    R.push_back(Y.size()-1);
    return R;
}

vector<double> inverse_function(vector<double> &X,vector<double> &Y,double y){
    vector<double> R(0);
    vector<double> K = inverse_solution_range(Y,y);
    vector<double> Z = inverse_monotonic_range(Y);
    int m = K.size()-1;
    int p = Z.size()-1;
    int f = 6;
    // cout << "K : \n";
    // for(auto i : K){
    //     cout << X[i] << ' ';
    // }
    // cout << '\n';

    // cout << "Z : \n";
    // for(auto i : Z){
    //     cout << i << ' ';
    // }
    // cout << '\n';

    int i = 0;
    int j;
    vector<double> U,V,L;
    double x, l, r;
    while(i <= m){
        // cout << i << '\n';
        j = 0;
        while(j < p){
            // cout << j << '\n';
            if (K[i] >= Z[j] && K[i]+1 <= Z[j+1])
            {
                if(Z[j+1] - Z[j] >= 5){
                    U.resize(0);
                    for(int q = Z[j];q <= Z[j+1];q++){
                        U.push_back(X[q]);
                    }
                    vector<int> M = newton_extract_points(U, X[K[i]], f);
                    l = M[0]+Z[j];
                    r = M[M.size()-1]+Z[j];
                }else{
                    l = Z[j];
                    r = Z[j+1];
                }

                // cout << X[l] << '-' << X[r] << '\n';

                U.resize(0);
                V.resize(0);
                for(int q = l;q <= r;q++){
                    U.push_back(Y[q]);
                    V.push_back(X[q]);
                }
                L = lagrange_non_equidistant(U,V);
                    // for(auto i : L){
                    //     cout << i << ' ';
                    // }
                    // cout << '\n';
                x = horner_polynominal_value(L,y);
                R.push_back(x);
            }
            j++;
        }
        i++;
    }
    return R; 
}


vector<double> inverse_newton_forward(vector<double> &X,vector<double> &Y,double y,double error){
    vector<double> K = inverse_solution_range(Y,y);
    int k = 11;
    cout << "K : \n";
    for(auto i : K){
        cout << i << ' ';
    }
    cout << '\n';

    vector<double> R(0);
    vector<double> U,V;
    double x, t_0, t_1;
    for(int i = 0;i < K.size();i++){
        // if (i == 0){
        //     continue;
        // }
        // cout << K[i] << "-\n";

        vector<int> M = newton_extract_points(X, X[K[i]], k);
        vector<double> X1(k);
        vector<double> Y1(k);

        for(auto i : M){
            cout << i << ' ';
        }
        cout << '\n';

        for(int i = 0;i < M.size();i++){
            X1[i] = X[M[i]];
            Y1[i] = Y[M[i]];
        }

        vector<vector<double>>  Delta = inverse_forward_diff(Y1);

        // cout << "Bang ti hieu : \n";
        // for(auto a : Delta){
        //     for(auto b : a){
        //         cout << b << ' ';
        //     }
        //     cout << '\n';
        // }

        int n = Y1.size()-1;

        U.resize(n+1);
        for(auto &i : U){
            i = 0;
        }
        U[0] = y - Y1[0];

        cout << "U : \n";
        for(auto i : U){
            cout << i << ' ';
        }
        cout << '\n';

        for(int q = 2;q <= n;q++){
            V = {1};
            for(int w = 0;w < q;w++){
                x = w;
                V = horner_polynominal_multiply(V,x);
            }

            cout << "V : \n";
            for(auto i : V){
                cout << i << ' ';
            }
            cout << '\n';

            for(int w = 0;w < V.size();w++){
                U[w] -= (V[w]*Delta[q-1][0])/factorial(q);
            }

            cout << "U : \n";
            for(auto i : U){
                cout << i << ' ';
            }
            cout << '\n';
        }

        for(int q = 0;q < U.size();q++){
            U[q] /= Delta[0][0];
        }
        t_1 = (y - Y1[0])/Delta[0][0];  \

        // for (int q = 0;q < 4;q++){
        //     t_0 = t_1;
        //     t_1 = horner_polynominal_value(U,t_0);
        //     cout << "t : \n" << t_0 << '-' << t_1 << '\n';
        // }

        do{
            t_0 = t_1;
            t_1 = horner_polynominal_value(U,t_0);

            cout << "t : \n" << t_0 << '-' << t_1 << '\n';
        }while(abs(t_0 - t_1) < error);
        cout << X1[0] << '-' << t_1 << '\n';
        x = X1[0] + t_1*(X1[1] - X1[0]);
        R.push_back(x);
        // cout << "======\n" << x << '\n';
    }
    return R; 
}

vector<double> inverse_newton_back(vector<double> &X,vector<double> &Y,double y,double error){
    vector<double> K = inverse_solution_range(Y,y);
    int k = 7;
    // cout << "K : \n";
    // for(auto i : K){
    //     cout << i << ' ';
    // }
    // cout << '\n';

    vector<double> R(0);
    vector<double> U,V;
    double x, t_0, t_1;
    for(int i = 0;i < K.size();i++){
        // if (i > 0){
        //     continue;
        // }
        // cout << K[i] << "-\n";

        vector<int> M = newton_extract_points(X, X[K[i]], k);
        vector<double> X1(k);
        vector<double> Y1(k);

        // for(auto i : M){
        //     cout << i << ' ';
        // }
        // cout << '\n';

        for(int i = 0;i < M.size();i++){
            X1[i] = X[M[i]];
            Y1[i] = Y[M[i]];
        }

        // for(auto i : Y1){
        //     cout << i << ' ';
        // }
        // cout << '\n';

        vector<vector<double>>  Delta = inverse_back_diff(Y1);

        // cout << "Bang ti hieu : \n";
        // for(auto a : Delta){
        //     for(auto b : a){
        //         cout << b << ' ';
        //     }
        //     cout << '\n';
        // }

        int n = Y1.size()-1;

        U.resize(n+1);
        for(auto &i : U){
            i = 0;
        }
        U[0] = y - Y1[n];

        // cout << "U : \n";
        // for(auto i : U){
        //     cout << i << ' ';
        // }
        // cout << '\n';

        for(int q = 2;q <= n;q++){
            V = {1};
            for(int w = 0;w < q;w++){
                x = -w; 
                V = horner_polynominal_multiply(V,x);
            }

            // cout << "V : \n";
            // for(auto i : V){
            //     cout << i << ' ';
            // }
            // cout << '\n';

            for(int w = 0;w < V.size();w++){
                U[w] -= (V[w]*Delta[q-1][n-q])/factorial(q);
            }

            // cout << "U : \n";
            // for(auto i : U){
            //     cout << i << ' ';
            // }
            // cout << '\n';
        }

        for(int q = 0;q < U.size();q++){
            U[q] /= Delta[0][n-1];
        }
        t_1 = (y - Y1[0])/Delta[0][n-1];  

        // for (int q = 0;q < 4;q++){
        //     t_0 = t_1;
        //     t_1 = horner_polynominal_value(U,t_0);
        //     cout << "t : \n" << t_0 << '-' << t_1 << '\n';
        // }

        do{
            t_0 = t_1;
            t_1 = horner_polynominal_value(U,t_0);
            cout << "t : \n" << t_0 << '-' << t_1 << '\n';
        }while(abs(t_0 - t_1) > error);
        // cout << X1[0] << '-' << t_1 << '\n';
        x = X1[n] + t_1*(X1[1] - X1[0]);
        R.push_back(x);
        // cout << "======\n" << x << '\n';
    }

    return R; 
}

// vector<double> inverse_newton_back(vector<double> &X,vector<double> &Y,double y,double error){
//     int n = X.size()-1;
//     vector<double> R(0);
//     vector<double> K = inverse_solution_range(Y,y);
//     vector<double> Z = inverse_monotonic_range(Y);
//     vector<vector<double>> Delta = inverse_back_diff(Y);
//     int m = K.size()-1;
//     int p = Z.size()-1;
//     int i = 0;
//     int j;
//     vector<double> U,V;
//     double x,t_0,t_1, s,f;
//     while(i <= m){
//         j = 0;
//         while(j < p){
//             if (K[i] >= Z[j] && K[i]+1 <= Z[j+1])
//             {
//                 s = Z[j];
//                 f = Z[j+1] - Z[j];
//                 U.resize(f+1,0);
//                 U[0] = y - Y[n];
//                 for(int q = 2;q <= n;q++){
//                     V = {1};
//                     for(int w = 0;w < q;w++){
//                         x = -w; 
//                         V = horner_polynominal_multiply(V,x);
//                     }
//                     for(int w = 0;w < V.size();w++){
//                         U[w] -= (V[w]*Delta[q-1][n-q])/factorial(q);
//                     }
//                 }
//                 for(int q = 0;q < U.size();q++){
//                     U[q] /= Delta[0][0];
//                 }
//                 t_1 = (y - Y[n])/Delta[0][n-1];  
//                 do{
//                     t_0 = t_1;
//                     for(int q = 0;q < U.size();q++){
//                         t_1 += pow(t_0,q)*U[q];
//                     }
//                 }while(abs(t_0 - t_1) < error);
//                 x = X[n] + t_1*(X[1] - X[0]);
//                 R.push_back(x);
//             }
//             j++;
//         }
//         i++;
//     }
//     return R; 
// }

vector<vector<double>> inverse_forward_diff(vector<double> &Y){
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

vector<vector<double>> inverse_back_diff(vector<double> &Y){
    int n = Y.size()-1;
    int i = 1;
    vector<double> A = Y;
    vector<double> B;
    vector<vector<double>> ans;
    while(i <= n){
        B.resize(n-i+1);
        for(int j = n;j >= i;j--){
            B[j-i] = (A[j-i+1] - A[j-i]);
        }
        ans.push_back(B);
        A = B;
        B.resize(0);
        i++;
    }
    return ans;
}