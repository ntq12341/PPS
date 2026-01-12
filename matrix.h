/*gồm các thao tác với ma trận 
-matrix(n) : khởi tạo ma trận đơn vị cấp n
-matrix(m,n) : khởi tạo ma trận cỡ mxn
-matrix(double arr[m][n]) : khởi tạo ma trận theo arr
-mat.inp(): nhập ma trận 
-mat.out(): in ma trận 
-mat.row : số hàng
-mat.column : số cột
-mat.mat : mảng nền của ma trận
-mat.t() : ma trận chuyển vị
-mat.det() : định thức ma trận
-mat.ij(i,j) : ma trận con ij
-mat.one() : chuẩn 1 của ma trận
-mat.inf() : chuẩn inf của ma trận
-mat.row_dom() : chéo trội hàng
-mat.column_dom() : chéo trội cột
-mat.swap_row(r1,r2) : hoán đổi 2 hàng r1 và r2
-mat.trape() : ma trận dạng bậc thang
-mat.rank() : hạng ma trận
-mat.inverse() : ma trận nghịch đảo
*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

bool search_matrixlib(vector<int> &v, int x) {
	for (int i : v) {
		if (i == x) {
			return true;
		}
	}
	return false;
}

class matrix {
public:
	int row, column;
	vector<vector<double>> mat;

	matrix(){
		row = 0;
		column = 0;
	}

	matrix(int n){
		this->row = n;
		this->column = n;
		mat.resize(n);
		for(int i = 0;i < n;i++){
			mat[i].resize(n);
			for(int j = 0;j < n;j++){
				if(j == i){
					mat[i][j] = 1;
				}else{
					mat[i][j] = 0;
				}
			}
		}
	}

	matrix(int row, int column) {
		this->row = row;
		this->column = column;
		mat.resize(row);
		for (auto& vec : mat) {
			vec.resize(column);
			for(auto &d : vec){
				d = 0;
			}
		}
	}

	matrix(vector<vector<double>> &v){
		this->row = v.size();
		this->column = v[0].size();
		this->mat = v;
	}

	void inp() {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				cin >> mat[i][j];
			}
		}
	}

	void out() {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				if (abs(mat[i][j]) < 1e-4) {
					cout << 0 << ' ';
				}
				else {
					cout << mat[i][j] << ' ';
				}
			}
			cout << endl;
		}
		cout << endl;
	}

	void operator=(const matrix& mat_2) {
		row = mat_2.row;
		column = mat_2.column;
		mat = mat_2.mat;
	}

	matrix operator+(const matrix& mat_2) {
		if (mat_2.row != this->row || mat_2.column != this->column) {
			cout << "Please check number of rows or columns again" << endl;
			exit(1);
		}
		else {
			matrix res(row, column);
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < column; j++) {
					res.mat[i][j] = mat_2.mat[i][j] + mat[i][j];
				}
			}
			return res;
		}
	}

	bool operator==(const matrix& mat_2){
		for(int i = 0;i < row;i++){
			for(int j = 0;j < column;j++){
				if(mat[i][j] != mat_2.mat[i][j]) return false;
			}
		}
		return true;
	}

	matrix operator-(const matrix& mat_2) {
		if (mat_2.row != this->row || mat_2.column != this->column) {
			cout << "Please check number of rows or columns again" << endl;
			exit(1);
		}
		else {
			matrix res(row, column);
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < column; j++) {
					res.mat[i][j] = -mat_2.mat[i][j] + mat[i][j];
				}
			}
			return res;
		}
	}

	matrix operator*(const double& n) {
		matrix res(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				res.mat[i][j] = n * mat[i][j];
			}
		}
		return res;
	}

	matrix operator*(const matrix& mat_2) {
		if (this->column != mat_2.row) {
			cout << "Please check number of rows or columns again" << endl;
			exit(1);
		}
		else {
			matrix res(this->row, mat_2.column);
			for (int i = 0; i < this->row; i++) {
				for (int j = 0; j < mat_2.column; j++) {
					for (int k = 0; k < this->column; k++) {
						res.mat[i][j] += (this->mat[i][k] * mat_2.mat[k][j]);
					}
				}
			}
			return res;
		}
	}

	matrix t() {
		matrix res(column, row);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				res.mat[j][i] = mat[i][j];
			}
		}
		return res;
	}

	matrix ij(int n, int m) {
		vector<vector<double>> ans(row - 1,vector<double> (column-1));
		for (int i = 0; i < row; i++) {
			if (i < n) {
				for (int j = 0; j < column; j++) {
					if (j < m) {
						ans[i][j] = mat[i][j];
					}	
					else if (j > m) {
						ans[i][j - 1] = mat[i][j];
					}
				}
			}
			else if (i > n) {
				for (int j = 0; j < column; j++) {
					if (j < m) {
						ans[i-1][j] = mat[i][j];
					}
					else if (j > m) {
						ans[i-1][j - 1] = mat[i][j];
					}	
				}
			}
		}
		matrix a(ans);
		return a;
    }

    double det() {
	    if (row != column) {
		    return NAN;
	    }
	    else {
			if (row == 2) {
				return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
			}	
			else if (row == 1) {
				return mat[0][0];
			}
			else {
				double ans = 0;
				for (int i = 0; i < row; i++) {
					if (i % 2 == 0) {
						ans += mat[0][i] * (this->ij(0,i)).det();
					}
					else {
						ans -= mat[0][i] * (this->ij(0,i)).det();
					}
				}
				return ans;
			}
	    }
    }

    double one() {
	double ans = 0;
	for (int i = 0; i < column; i++) {
		double k = 0;
		for (int j = 0; j < row; j++) {
			k += abs(mat[j][i]);
		}
		ans = (k > ans) ? k : ans;
	}
	return ans;
    }

    double inf() {
	double ans = 0;
	for (int i = 0; i < row; i++) {
		double k = 0;
		for (int j = 0; j < column; j++) {
			k += abs(mat[i][j]);
		}
		ans = (k > ans) ? k : ans;
	}
	return ans;
    }

    bool row_dom(){
    	for(int i = 0;i < row;i++){
    		double sum = 0;
    		for(int j = 0;j < column;j++){
    			sum += abs(mat[i][j]);
    		}
    		if(sum >= 2*abs(mat[i][i])) return false;
    	}
    	return true;
    }

    bool column_dom(){
    	for(int j = 0;j < column;j++){
    		double sum = 0;
    		for(int i = 0;i < row;i++){
    			sum += abs(mat[i][j]);
    		}
    		if(sum >= 2*abs(mat[j][j])) return false;
    	}
    	return true;
    }

    void swap_row(int r1,int r2){
		for (int i = 0; i < column; i++) {
			swap(mat[r1][i], mat[r2][i]);
		}
    }

    matrix trape() {
    	vector<vector<double>> v = this->mat;
		int r = min(row,column);
		int i = 0, j = 0;
		while (i < r && j < r)
		{
			if (v[i][j] != 0) {
			for (int k = i + 1; k < r; k++) {
				double x = v[k][j] / v[i][j];
				for (int h = 0; h < r; h++) {
					v[k][h] -= x * v[i][h];
				}
			}
			i++;
			j++;
			}
			else {
				if (i == r - 1) {
					break;
				}
				int t = i + 1;
				while (v[t][j] == 0) {
					t++;
					if (t == r ) {
						break;
					}
				}
				if (t < r) {
					for (int z = 0;z < column;z++) {
						swap(v[t][z], v[i][z]);
					}
				}
				else {
					j++;
					if (j == r) {
						break;
					}
				}
			}
		}
		matrix ans(v);
		return ans;
    }

    int rank() {
		vector<int> zero(row);
		(this->trape()).out();
		vector<vector<double>> v = (this->trape()).mat;
		for (int i = 0; i < row; i++) {
			zero[i] = -1;
			for (int j = 0; j < column; j++) {
				if (v[i][j] != 0) {
					zero[i] = j;
					break;
				}
			}
		}
		int rank = row;
		int k = row - 1;
		while ( k >= 0 && zero[k] == -1) {
			k--;
			rank--;
		}
		return rank;
    }

    double norm_euclide(){
    	double ans = 0;
		for (int i = 0; i < column; i++) {
			for (int j = 0; j < row; j++) {
				ans += mat[i][j]*mat[i][j];
			}	
		}
		return sqrt(ans);
    }

    // matrix inverse(){
    // matrix b(row);
    // matrix a(row,column);
    // a.mat = this->mat;
    // vector<vector<int>> lis(2);
    // int i,j;
    // bool check;
	// bool check1;
	// double max;
	// do {
	// 	max = 0;
	// 	check = false;
	// 	check1 = false;
	// 	for (int q = 0; q < row; q++) {
	// 		if (!search_matrixlib(lis[0], q)) {
	// 			for (int w = 0; w < row; w++) {
	// 				if (!search_matrixlib(lis[1], w)) {
	// 					if (abs(a.mat[q][w]) == 1) {
	// 						i = q;
	// 						j = w;
	// 						lis[0].push_back(i);
	// 						lis[1].push_back(j);
	// 						check = true;
	// 						check1 = true;
	// 						break;
	// 					}
	// 					else {
	// 						if (abs(a.mat[q][w]) > max) {
	// 							max = abs(a.mat[q][w]);
	// 							i = q;
	// 							j = w;
	// 							check = true;
	// 						}
	// 					}
	// 				}
	// 			}
	// 			if (check1) {
	// 				break;
	// 			}
	// 		}
	// 	}
	// 	if (max > 0 && !check1) {
	// 		lis[0].push_back(i);
	// 		lis[1].push_back(j);
	// 	}
	// 	if (check) {
	// 		for (int q = 0; q < row; q++) {
	// 			if (q != i) {
	// 				double x = a.mat[q][j] / a.mat[i][j];
	// 				for (int w = 0; w < row; w++) {
	// 					a.mat[q][w] -= x * a.mat[i][w];
	// 				}
	// 				for (int w = 0; w < row; w++) {
	// 					b.mat[q][w] -= x * b.mat[i][w];
	// 				}
	// 			}
	// 			else {
	// 				double x = a.mat[i][j];
	// 				for (int w = 0; w < row; w++) {
	// 					a.mat[q][w] /= x;
	// 				}
	// 				for (int w = 0; w < row; w++) {
	// 					b.mat[q][w] /= x;
	// 				}
	// 			}
	// 		}
	// 	}
	// } while (check);
	// vector<int> ind(row,-1);
	// for (int k = 0; k < row; k++) {
	// 	for (int h = 0; h < 2*row; h++) {
	// 		if (h < row) {
	// 			if (a.mat[k][h] != 0) {
	// 				ind[k] = h;
	// 				break;
	// 			}
	// 		}
	// 		else {
	// 			if (b.mat[k][h - row] != 0) {
	// 				ind[k] = h;
	// 				break;
	// 			}
	// 		}
	// 	}
	// }
	// for(int p = 0;p < row - 2;p++){
	// 	bool sort = true;
	// 	for(int q = 0;q < row - p - 1;q++){
	// 		if(ind[q] > ind[q+1]){
	// 			a.swap_row(q,q+1);
	// 			b.swap_row(q,q+1);
	// 			swap(ind[q],ind[q+1]);
	// 			sort = false;
	// 		}
	// 	}
	// 	if(sort) break;
	// }
	// return b;
    // }

    matrix inverse(){
		int n = this->row;
		matrix a(n,n);
		a.mat = this->mat;
		// m = a^t * a 
		matrix m(n,n);
		m = a.t()*a;
		int i = 1;
		matrix m1;
		while(i <= n){
			matrix m2(i,i);
			if(i == 1){
				m2.mat[0][0] = 1/m.mat[0][0];
			}else{
				matrix a1(i-1,1); // ma trận cột ngoài của A
				matrix a2(1,i-1); // ma trận hàng dưới của A 
				for(int j = 0;j < i-1;j++){
					a1.mat[j][0] = m.mat[j][i-1];
				}
				a2 = a1.t();
				double a3 = m.mat[i-1][i-1]; // phần tử dưới cùng bên phải của A

				double b3 = 1/(a3 - ((a2*m1*a1).mat[0][0])); // phần tử dưới cùng bên phải của A^-1
				matrix b0(i-1,i-1); // ma trận trên cùng bên trái của A^-1
				matrix b1(i-1,1); // ma trận cột ngoài của A^-1
				matrix b2(1,i-1); // ma trận hàng dưới của A^-1
				matrix e(i-1); // ma trận đơn vị cấp i-1
				// b1 = -(m1*a1)/b3;
				// b2 = -(a2*m1)/b3;
				// b0 = m*(e + (a1*a2*m1)/b3);

				b0 = m1*(e + (a1*a2*m1)*b3);
				b1 = (m1*a1)*(-b3);
				b2 = (a2*m1)*(-b3);
				for(int j = 0;j < i-1;j++){
					m2.mat[j][i-1] = b1.mat[j][0];
					m2.mat[i-1][j] = b2.mat[0][j];
					for(int k = 0;k < i-1;k++){
						m2.mat[j][k] = b0.mat[j][k];
					}
				}
				m2.mat[i-1][i-1] = b3;
			}
			m1 = m2;
			i++;
		}
		matrix ans(n,n);
		ans = m1*a.t();
		return ans; 
	}
};
