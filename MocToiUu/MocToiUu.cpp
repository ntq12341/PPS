#include <iostream>
#include <vector>
#include<cmath>
#define PI 3.14159265359
using namespace std;

typedef long long ll;

//------------------------------------------------------------------------


//------------------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(){
    int n;
    double a,b;
    cin >>n>> a>> b;
    double z;
    vector<double> X(n+1,0);
    for(int i = 0;i <= n;i++){
        z = cos(((2*i+1)/(double)(2*n))*PI);
        // X[i] = (2*z - b-a)/(b-a);
        X[i] = (z*(b-a) + b+a)/2;
    }
    for(auto i : X){
        cout << i << '\n';
    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++