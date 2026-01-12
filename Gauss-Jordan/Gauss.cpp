#include "matrix.h"

// Quy trình thuận pp Gauss
void Gauss1(matrix &A,matrix &b){
	int m = A.row; // số hàng ma trận
	int nA = A.column; // số cột ma trận A
	int nb = b.column;// số cột ma trận b
	int l = (m <= nA) ? m : nA; // số lần lặp 
	for(int i = 0; i < l;i++){
		// TH1 : A_ii khác 0
		if(A.mat[i][i] != 0){
			for(int j = i+1; j < m;j++){
				double q = A.mat[j][i]/A.mat[i][i];
				for(int k = 0;k < nA;k++){
					A.mat[j][k] -= q * A.mat[i][k];
				}
				for(int k = 0;k < nb;k++){
					b.mat[j][k] -= q * b.mat[i][k];
				}
			}
		}
		// TH2 : A_ii = 0
		else{
			int f = i+1;
			bool check = true;
			while(A.mat[f][f] == 0) {
			    f++;
			    if(f == m){
			    	check = false;
			    	break;
			    }
			}
			if(check){
				A.swap_row(i,f);
				b.swap_row(i,f);
			}
		}
	}
}
// Quy trình nghịch pp Gauss trường hợp không suy biến
void Gauss2(matrix &A,matrix &b){
	int m = A.row; // số hàng ma trận
	int nA = A.column; // số cột ma trận A
	int nb = b.column;// số cột ma trận b
	for(int i = m-1;i >= 0;i--){
		for(int j = 0;j < nb;j++){
			for(int k = m-1;k >= i+1;k--){
				b.mat[i][j] -= A.mat[i][k]*b.mat[k][j];
			}
			b.mat[i][j] /= A.mat[i][i];
		}
	}
}

int main(){
	int m,nA,nb;
	cin >> m >> nA >> nb;
	matrix A(m,nA);
	matrix b(m,nb);
	A.inp();
	b.inp();
	Gauss1(A,b);
	A.out();
	b.out();
	Gauss2(A,b);
	A.out();
	b.out();
}