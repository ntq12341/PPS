// #include "Gauss-Jordan/Gauss-Jordan.cpp"
// #include "Lagrange/Lagrange.cpp" 
// #include "DaoHam/DaoHam.cpp" 
<<<<<<< HEAD
// #include "Euler/Euler2.cpp"
// #include "RK\RK.cpp"
// #include "Adam\Adam.cpp"  
#include "SaiPhan\SaiPhan.cpp" 
=======
#include "Euler/Euler2.cpp" 
>>>>>>> c4205160f93fe26a35c58aa1c58fa58ad919a256
#include <iomanip> 

#define PI 3.14159265359
#define E 2.718281828

int main(){
<<<<<<< HEAD
    double x0 = 0, xn = 3, a = 1, b = 0.25, h = 0.1, o1 = 2, o2 = 0.1;
    vector<double> SP = SaiPhan3(x0, xn, h, a, b, o1, o2); 
    for(auto u : SP){
        cout << u << '\n';
    }
    // for(int i = 0;i < Adam[0].size();i++){
    //     cout << Adam[0][i] << " : " << Adam[1][i] << '\n';
    // }
=======
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
>>>>>>> c4205160f93fe26a35c58aa1c58fa58ad919a256
}
