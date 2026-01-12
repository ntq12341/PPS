// #include "Gauss-Jordan/Gauss-Jordan.cpp"
// #include "Lagrange/Lagrange.cpp" 
// #include "DaoHam/DaoHam.cpp" 
#include "Euler/Euler2.cpp" 
#include <iomanip> 

#define PI 3.14159265359
#define E 2.718281828

int main(){
    cout << setprecision(10) ;
    double x0 = 12;
    double y0 = 18;
    double z0 = 8;
    double t0 = 0;
    double tn = 1500;
    double h = 0.1; 

    vector<vector<double>> A = Euler_Hien_Hept(x0,y0,z0,t0,tn,h);
    for(int i = 0;i < A[0].size();i++){
        cout << A[0][i] << ':' << A[1][i] << ':' << A[2][i] << ':' << A[3][i] << '\n';
    }
}
