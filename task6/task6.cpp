#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Eigen/Dense"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


void print_M(double **matrix, int M, int N){
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
			cout << setw(13) << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


void print_V(double *vec, int M){
	for (int i = 0; i < M; i++){
		cout << setw(13) << vec[i] << endl;
	}
	cout << endl;
}


double f(double x){
	return sin(x) * x;
	//return x * x * x;
}


double Phi_i(double* denoms, int size, double x, double a, double step, int i, int k){
	double value_x = 0.0, num = 1.0;
	for (int j = 0; j < size; j++){
		if (j != (i - k * (size - 1))){
			num *= (x - (a + (j + k * (size - 1)) * step));
		}
	}
	value_x = num / denoms[i - k * (size - 1)];
	return value_x;
}


double* Gauss(double **A, double *B, int M){
	double    *x = new double[M];
	double  *B_t = new double[M];
	double **A_t = new double*[M];
	double tmp = 0.;
	
	for (int i = 0; i < M; i++){
		A_t[i] = new double[M];
		B_t[i] = B[i];
		for (int j = 0; j < M; j++){
			A_t[i][j] = A[i][j];
		}
	}
	
	for (int i = 0; i < M; i++){
		for (int j = i + 1; j < M; j++){
			tmp = A_t[j][i]/A_t[i][i];
			
			for (int k = 0; k < M; k++){
				A_t[j][k] -= A_t[i][k] * tmp;
			}
			B_t[j] -= tmp * B_t[i];
		}
		
	}
	
	tmp = 0.;
	for (int i = M - 1 ; i >= 0; i--){
		for (int j = M - 1; j > i; j--){
			tmp += A_t[i][j] * x[j];
		}
		x[i] = (B_t[i] - tmp)/A_t[i][i];
		tmp = 0.;
	}
	
	//print_M(A, M, M);
	
	for (int i = 0; i < M; i++){
		delete [] A_t[i];
	}
	delete [] A_t;
	delete [] B_t;
	
	
	return x;
}


double* LU(double **A, double *B, int M){
	double  *x = new double[M];
	double  *y = new double[M];
	double **L = new double*[M];
	double **U = new double*[M];
	
	
	for (int i = 0; i < M; i++){
		L[i] = new double[M];
		U[i] = new double[M];
	}
	
	
	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++){
			U[i][j] = 0.;
			L[i][j] = 0.;
		}
		L[i][i] = 1.;
	}
	
	
	double tmp_1 = 0., tmp_2 = 0.;
	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++){
			if (i <= j){
				for (int k = 0; k < i; k++){
					tmp_1 += L[i][k] * U[k][j];
				}
				U[i][j] = A[i][j] - tmp_1;
				tmp_1 = 0.;
			}
			else if (i > j){
				for (int k = 0; k < j; k++){
					tmp_2 += L[i][k] * U[k][j];
				}
				L[i][j] = (A[i][j] - tmp_2)/U[j][j];
				tmp_2 = 0.;
			}
		}
	}
	
	double tmp = 0.;
	for (int i = 0 ; i < M; i++){
		for (int j = 0; j < i; j++){
			tmp += L[i][j] * y[j];
		}
		y[i] = B[i] - tmp;
		tmp = 0.;
	}
	
	
	for (int i = M - 1 ; i >= 0; i--){
		for (int j = M - 1; j > i; j--){
			tmp += U[i][j] * x[j];
		}
		x[i] = (y[i] - tmp)/U[i][i];
		tmp = 0.;
	}
	
	//print_M(L, M, M);
	//print_M(U, M, M);
	
	
	for (int i = 0; i < M; i++){
		delete [] L[i];
		delete [] U[i];
	}
	delete [] L;
	delete [] U;
	delete [] y;
	
	return x;
}


double* Cholesky(double **A, double *B, int M){
	double  *x = new double[M];
	double  *y = new double[M];
	double **L = new double*[M];
	
	
	for (int i = 0; i < M; i++){
		L[i] = new double[M];
	}
	
	
	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++){
			L[i][j] = 0.;
		}
	}
	
	
	double tmp = 0.;
	for (int i = 0; i < M; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k < j; k++){
				tmp += L[i][k] * L[j][k];
			}
			
			if (i == j){
				L[i][j] = sqrt(A[i][i] - tmp);
			}
			else{
				L[i][j] = (A[i][j] - tmp)/L[j][j];
			}
			tmp = 0.;
		}
	}
	
	//print_M(L, M, M);
	
	//tmp = 0.;
	for (int i = 0 ; i < M; i++){
		for (int j = 0; j < i; j++){
			tmp += L[i][j] * y[j];
		}
		y[i] = (B[i] - tmp)/L[i][i];
		tmp = 0.;
	}
	
	
	for (int i = M - 1 ; i >= 0; i--){
		for (int j = M - 1; j > i; j--){
			tmp += L[j][i] * x[j];
		}
		x[i] = (y[i] - tmp)/L[i][i];
		tmp = 0.;
	}
	
	
	for (int i = 0; i < M; i++){
		delete [] L[i];
	}
	delete [] L;
	delete [] y;
	
	
	return x;
}


double* Relax(double **A, double *B, int M){
	double* x = new double[M];
	double omega = 1.2, eps = 1e-5, r = 1.;
	int n = 0;
	
	
	for (int i = 0; i < M; i++){
		x[i] = 0.;
	}
	
	
	double tmp = 0.;
	while (r > eps and n < 50){
		r = 0.;
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				tmp += A[i][j] * x[j];
			}
			x[i] += omega * (B[i] - tmp)/A[i][i];
			
			r += (B[i] - tmp) * (B[i] - tmp);
				
			tmp = 0.;
		}
		r = sqrt(r);
		n++;
	}
	
	cout << endl << "Relax    |  Number of iterations: " << n << endl;
	
	return x;
}


double* Gradient(double **A, double *B, int M){
	double* x  = new double[M];
	double* r  = new double[M];
	double* p  = new double[M];
	double* Ap = new double[M];
	double eps = 1e-5, r_norm = 1., alpha, beta;
	int n = 0;
	
	
	for (int i = 0; i < M; i++){
		x[i] = 0.;
		r[i] = B[i];
		p[i] = r[i];
	}
	
	
	double tmp_1 = 0., tmp_2 = 0., tmp_3 = 0., tmp = 0.;						// tmp_1 = (r_k, r_k), tmp_2 = (Ap, p), tmp_3 = (r_k+1, r_k+1)
	while (r_norm > eps and n <= M){
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				tmp += A[i][j] * p[j];
			}
			
			Ap[i] = tmp;
			tmp = 0.;
			
			tmp_1 += r[i] * r[i];
			tmp_2 += p[i] * Ap[i];
		}
			
		alpha = tmp_1/tmp_2;
		
		for (int i = 0; i < M; i++){	
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			
			tmp_3 += r[i] * r[i];
		}
		
		beta = tmp_3/tmp_1;
		
		for (int i = 0; i < M; i++){			
			p[i] = beta * p[i] + r[i];
		}
		
		r_norm = tmp_3;
		r_norm = sqrt(r_norm);
		
		tmp_1 = 0.;
		tmp_2 = 0.;
		tmp_3 = 0.;
		
		n++;
	}
	
	//cout << r_norm << endl;
	
	cout << "Gradient |  Number of iterations: " << n << endl;
	
	delete [] p;
	delete [] r;
	delete [] Ap;
	
	return x;
}


int main(void){
	double a = -3.14, b = 3.14, h, len, h_total, den = 1.0, delta;
	double*  points;
	double*  denoms;
	double** A_b;
	double*  B_b;
	double*  X_gs;
	double*  X_lu;
	double*  X_ch;
	double*  X_rx;
	double*  X_gr;
	int K = 5, N = 4, M, L = N + 7, Z = 300;
	
	M       = K * (N - 1) + 1;
	h       = (b - a) / (M - 1);
	h_total = pow(h, N - 1);
	len     = (b - a) / K;	
	delta   = (b - a) / (Z - 1);
	
	MatrixXd A = MatrixXd::Zero(M, M);
	VectorXd B = VectorXd::Zero(M);
	points = new double[L * K];
	denoms = new double[N];
	
	
	uniform_real_distribution<double> distribution{a, a + len};
	mt19937 generator{random_device().operator ()()};
	
	cout.setf(ios::scientific);
	

	for (int k = 0; k < K; k++){
		for (int i = 0; i < L; i++){
			points[i + k * L] = len * k + distribution(generator);
		}
	}

	sort(points, points + K * L);
	
	for (int j = 0; j < N; j++){																					// Заполнение массива знаменателей
		for (int k = 0; k < N; k++){
			if (k != j)
				den *= (j - k);
		}
		denoms[j] = den * h_total;
		den = 1.0;
	}
	
	A_b = new double*[M];
	B_b = new double[M];
	
	for (int i = 0; i < M; i++){
		B_b[i] = 0.;
		A_b[i] = new double[M];
		for (int j = 0; j < M; j++){
			A_b[i][j] = 0.;
		}
	}
	
	double tmp1 = 0.0, tmp2 = 0.0;																			// Заполнение матрицы A и вектора B
	for (int k = 0; k < K; k++){
		for (int i = k * (N - 1); i <= (k + 1) * (N - 1); i++){
			for (int j = i; j <= (k + 1) * (N - 1); j++){
				tmp1 = 0.0;
				for (int n = k * L; n < (k + 1) * L; n++){
					tmp1 += Phi_i(denoms, N, points[n], a, h, j, k) * Phi_i(denoms, N, points[n], a, h, i, k);
				}
				A_b[i][j] += tmp1;
				A(i, j) += tmp1;
				if (i != j){
					A(j, i) += tmp1;
					A_b[j][i] += tmp1;
				}
			}
			for (int m = k * L; m < (k + 1) * L; m++){
				tmp2 += Phi_i(denoms, N, points[m], a, h, i, k) * f(points[m]);
			}
			B(i)   += tmp2;
			B_b[i] = B(i);			
			tmp2 = 0.0;
		}
	}
	
	//cout << endl << A << endl << endl << B << endl << endl;
	VectorXd C = A.colPivHouseholderQr().solve(B);										// Решение системы с помощью библиотечной функции
	
	
	
	/*
	int k = 0;
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
			k = i + j;
			if (k >= 0 and k < M){
				A_b[i][j] = A(i, k);
			}
			else
				A_b[i][j] = 0;			
			
		}
	}*/
	
	
	//cout << endl << C << endl << endl;
	
	X_lu = LU(A_b, B_b, M);
	//print_V(X_lu, M);
	
	
	X_ch = Cholesky(A_b, B_b, M);
	//print_V(X_ch, M);
	
	
	X_rx = Relax(A_b, B_b, M);
	//print_V(X_rx, M);
		
	
	X_gr = Gradient(A_b, B_b, M);
	//print_V(X_gr, M);
		
	
	X_gs = Gauss(A_b, B_b, M);
	//print_V(X_gs, M);
	
	
	double L_1_ch = 0., L_1_lu = 0., L_1_gs = 0., L_1_gr = 0., L_1_rx = 0.;
	double L_2_ch = 0., L_2_lu = 0., L_2_gs = 0., L_2_gr = 0., L_2_rx = 0.;
	double L_i_ch = 0., L_i_lu = 0., L_i_gs = 0., L_i_gr = 0., L_i_rx = 0.;
	double L_1_C  = 0., L_2_C  = 0., L_i_C  = 0.;
	double L_1_B  = 0., L_2_B  = 0., L_i_B  = 0.;
	double r_1_ch = 0., r_1_lu = 0., r_1_gs = 0., r_1_gr = 0., r_1_rx = 0.;
	double r_2_ch = 0., r_2_lu = 0., r_2_gs = 0., r_2_gr = 0., r_2_rx = 0.;
	double r_i_ch = 0., r_i_lu = 0., r_i_gs = 0., r_i_gr = 0., r_i_rx = 0.;
	double tmp_1  = 0., tmp_2  = 0., tmp_3  = 0., tmp_4  = 0., tmp_5  = 0.;
	
	for (int i = 0; i < M; i++){
		L_1_ch += abs(X_ch[i] - C(i));
		L_1_gr += abs(X_gr[i] - C(i));
		L_1_gs += abs(X_gs[i] - C(i));
		L_1_lu += abs(X_lu[i] - C(i));
		L_1_rx += abs(X_rx[i] - C(i));
		
		L_2_ch += (X_ch[i] - C(i)) * (X_ch[i] - C(i));
		L_2_gr += (X_gr[i] - C(i)) * (X_gr[i] - C(i));
		L_2_gs += (X_gs[i] - C(i)) * (X_gs[i] - C(i));
		L_2_lu += (X_lu[i] - C(i)) * (X_lu[i] - C(i));
		L_2_rx += (X_rx[i] - C(i)) * (X_rx[i] - C(i));
		
		L_i_ch = max(L_i_ch, abs(X_ch[i] - C(i)));
		L_i_gr = max(L_i_gr, abs(X_gr[i] - C(i)));
		L_i_gs = max(L_i_gs, abs(X_gs[i] - C(i)));
		L_i_lu = max(L_i_lu, abs(X_lu[i] - C(i)));
		L_i_rx = max(L_i_rx, abs(X_rx[i] - C(i)));
		
		L_1_C += abs(C(i));
		L_2_C += C(i) * C(i);
		L_i_C = max(L_i_C, abs(C(i)));
		
		L_1_B += abs(B_b[i]);
		L_2_B += B_b[i] * B_b[i];
		L_i_B = max(L_i_B, abs(B_b[i]));
		
		
		for (int j = 0; j < M; j++){
			tmp_1 += A_b[i][j] * X_ch[j];
			tmp_2 += A_b[i][j] * X_gr[j];
			tmp_3 += A_b[i][j] * X_gs[j];
			tmp_4 += A_b[i][j] * X_lu[j];
			tmp_5 += A_b[i][j] * X_rx[j];
		}
		
		r_1_ch += abs(tmp_1 - B_b[i]);
		r_1_gr += abs(tmp_2 - B_b[i]);
		r_1_gs += abs(tmp_3 - B_b[i]);
		r_1_lu += abs(tmp_4 - B_b[i]);
		r_1_rx += abs(tmp_5 - B_b[i]);
		
		r_2_ch += (tmp_1 - B_b[i]) * (tmp_1 - B_b[i]);
		r_2_gr += (tmp_2 - B_b[i]) * (tmp_2 - B_b[i]);
		r_2_gs += (tmp_3 - B_b[i]) * (tmp_3 - B_b[i]);
		r_2_lu += (tmp_4 - B_b[i]) * (tmp_4 - B_b[i]);
		r_2_rx += (tmp_5 - B_b[i]) * (tmp_5 - B_b[i]);
		
		r_i_ch = max(r_i_ch, abs(tmp_1 - B_b[i]));
		r_i_gr = max(r_i_gr, abs(tmp_2 - B_b[i]));
		r_i_gs = max(r_i_gs, abs(tmp_3 - B_b[i]));
		r_i_lu = max(r_i_lu, abs(tmp_4 - B_b[i]));
		r_i_rx = max(r_i_rx, abs(tmp_5 - B_b[i]));
		
		tmp_1 = 0.;
		tmp_2 = 0.;
		tmp_3 = 0.;
		tmp_4 = 0.;
		tmp_5 = 0.;
	}
	
	L_2_ch = sqrt(L_2_ch);
	L_2_lu = sqrt(L_2_lu);
	L_2_gs = sqrt(L_2_gs);
	L_2_gr = sqrt(L_2_gr);
	L_2_rx = sqrt(L_2_rx);
	
	r_2_ch = sqrt(r_2_ch);
	r_2_gr = sqrt(r_2_gr);
	r_2_gs = sqrt(r_2_gs);
	r_2_lu = sqrt(r_2_lu);
	r_2_rx = sqrt(r_2_rx);
	
	
	cout << "         |-----------------------------------------------------------------------------------------|" << endl;
	cout << "         |                  Absolute                  |                  Relative                  |" << endl;
	cout << "         |     L_1      |      L_2     |      L_i     |     L_1      |      L_2     |      L_i     |" << endl;
	cout << "Gauss    |" << setw(14) << L_1_gs << "|" << setw(14) << L_2_gs << "|" << setw(14) << L_i_gs << "|";
	cout << setw(14) << L_1_gs/L_1_C << "|" << setw(14) << L_2_gs/L_2_C << "|" << setw(14) << L_i_gs/L_i_C << "|" << endl;
	cout << "LU       |" << setw(14) << L_1_lu << "|" << setw(14) << L_2_lu << "|" << setw(14) << L_i_lu << "|";
	cout << setw(14) << L_1_lu/L_1_C << "|" << setw(14) << L_2_lu/L_2_C << "|" << setw(14) << L_i_lu/L_i_C << "|" << endl;
	cout << "Cholesky |" << setw(14) << L_1_ch << "|" << setw(14) << L_2_ch << "|" << setw(14) << L_i_ch << "|";
	cout << setw(14) << L_1_ch/L_1_C << "|" << setw(14) << L_2_ch/L_2_C << "|" << setw(14) << L_i_ch/L_i_C << "|" << endl;
	cout << "Relax    |" << setw(14) << L_1_rx << "|" << setw(14) << L_2_rx << "|" << setw(14) << L_i_rx << "|";
	cout << setw(14) << L_1_rx/L_1_C << "|" << setw(14) << L_2_rx/L_2_C << "|" << setw(14) << L_i_rx/L_i_C << "|" << endl;
	cout << "Gradient |" << setw(14) << L_1_gr << "|" << setw(14) << L_2_gr << "|" << setw(14) << L_i_gr << "|";
	cout << setw(14) << L_1_gr/L_1_C << "|" << setw(14) << L_2_gr/L_2_C << "|" << setw(14) << L_i_gr/L_i_C << "|" << endl;
	cout << "         |-----------------------------------------------------------------------------------------|" << endl;
	cout << "         |                                         Residue                                         |" << endl;
	cout << "         |                  Absolute                  |                  Relative                  |" << endl;
	cout << "         |     L_1      |      L_2     |      L_i     |     L_1      |      L_2     |      L_i     |" << endl;
	cout << "Gauss    |" << setw(14) << r_1_gs << "|" << setw(14) << r_2_gs << "|" << setw(14) << r_i_gs << "|";
	cout << setw(14) << r_1_gs/L_1_B << "|" << setw(14) << r_2_gs/L_2_B << "|" << setw(14) << r_i_gs/L_i_B << "|" << endl;
	cout << "LU       |" << setw(14) << r_1_lu << "|" << setw(14) << r_2_lu << "|" << setw(14) << r_i_lu << "|";
	cout << setw(14) << r_1_lu/L_1_B << "|" << setw(14) << r_2_lu/L_2_B << "|" << setw(14) << r_i_lu/L_i_B << "|" << endl;
	cout << "Cholesky |" << setw(14) << L_1_ch << "|" << setw(14) << L_2_ch << "|" << setw(14) << L_i_ch << "|";
	cout << setw(14) << r_1_ch/L_1_B << "|" << setw(14) << r_2_ch/L_2_B << "|" << setw(14) << r_i_ch/L_i_B << "|" << endl;
	cout << "Relax    |" << setw(14) << r_1_rx << "|" << setw(14) << r_2_rx << "|" << setw(14) << r_i_rx << "|";
	cout << setw(14) << r_1_rx/L_1_B << "|" << setw(14) << r_2_rx/L_2_B << "|" << setw(14) << r_i_rx/L_i_B << "|" << endl;
	cout << "Gradient |" << setw(14) << r_1_gr << "|" << setw(14) << r_2_gr << "|" << setw(14) << r_i_gr << "|";
	cout << setw(14) << r_1_gr/L_1_B << "|" << setw(14) << r_2_gr/L_2_B << "|" << setw(14) << r_i_gr/L_i_B << "|" << endl;
	cout << "         |-----------------------------------------------------------------------------------------|" << endl;
	
	
	
	
	
	for (int i = 0; i < M; i++)
		delete [] A_b[i];
	delete [] A_b;
	delete [] points;
	delete [] denoms;
	delete [] B_b;
	delete [] X_ch;
	delete [] X_lu;
	delete [] X_gr;
	delete [] X_gs;
	delete [] X_rx;
	
	
	return 0;
}