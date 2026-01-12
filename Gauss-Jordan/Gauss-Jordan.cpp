#include "../matrix.h"

matrix gauss_jordan(matrix &a_i, matrix &b_i){
    int row = a_i.row;
    //matrix b(a.row);
    matrix a = a_i;
    matrix b = b_i;
    vector<vector<int>> lis(2);
    int i,j;
    bool check;
	bool check1;
	double max;
	do {
		max = 0;
		check = false;
		check1 = false;
		for (int q = 0; q < row; q++) {
			if (!search_matrixlib(lis[0], q)) {
				for (int w = 0; w < row; w++) {
					if (!search_matrixlib(lis[1], w)) {
						// tìm |a_ij| = 1
						if (abs(a.mat[q][w]) == 1) {
							i = q;
							j = w;
							lis[0].push_back(i);
							lis[1].push_back(j);
							check = true;
							check1 = true;
							break;
						}
						// tìm max|a_ij|
						else {
							if (abs(a.mat[q][w]) > max) {
								max = abs(a.mat[q][w]);
								i = q;
								j = w;
								check = true;
							}
						}
					}
				}
				// tìm đc |a_ij| = 1 thì quyết định luôn phần tử khử
				if (check1) {
					break;
				}
			}
		}
		// loại hàng và cột chứa phần tử khử khỏi danh sách duyệt vòng sau
		if (max > 0 && !check1) {
			lis[0].push_back(i);
			lis[1].push_back(j);
		}
		if (check) {
			for (int q = 0; q < row; q++) {
				if (q != i) {
					double x = a.mat[q][j] / a.mat[i][j];
					for (int w = 0; w < a.column; w++) {
						a.mat[q][w] -= x * a.mat[i][w];
					}
					for (int w = 0; w < b.column; w++) {
						b.mat[q][w] -= x * b.mat[i][w];
					}
				}
				else {
					double x = a.mat[i][j];
					for (int w = 0; w < a.column; w++) {
						a.mat[q][w] /= x;
					}
					for (int w = 0; w < b.column; w++) {
						b.mat[q][w] /= x;
					}
				}
			}
		}
	} while (check);
	vector<int> ind(row,-1);
	// mảng ind ghi lại vị trí phần tử khác 0 đầu tiên mỗi hàng
	for (int k = 0; k < row; k++) {
		for (int h = 0; h < a.column + b.column; h++) {
			if (h < a.column) {
				if (a.mat[k][h] != 0) {
					ind[k] = h;
					break;
				}
			}
			else {
				if (b.mat[k][h - a.column] != 0) {
					ind[k] = h;
					break;
				}
			}
		}
	}
	// sắp xếp lại ind theo thứ tự tăng
	for(int p = 0;p < row - 1;p++){
		bool sort = true;
		for(int q = 0;q < row - p - 1;q++){
			if(ind[q] > ind[q+1]){
				a.swap_row(q,q+1);
				b.swap_row(q,q+1);
				swap(ind[q],ind[q+1]);
				sort = false;
			}
		}
		if(sort) break;
	}
	// a.out();
	// b.out();
	return b;
}

// int main(){
// 	matrix a(5);
// 	matrix b(5,1);
// 	a.inp();
// 	b.inp();
// 	matrix x = gauss_jordan(a,b);
// 	// a.out();
// 	// b.out();
// 	x.out(); 
// }