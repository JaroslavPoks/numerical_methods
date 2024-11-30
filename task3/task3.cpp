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


int main(void){
	double a = -3.14, b = 3.14, h, len, h_total, den = 1.0, delta;
	double L_1_f = 0.0, L_2_f = 0.0, L_inf_f = 0.0;
	double L_1_abs_r = 0.0, L_2_abs_r = 0.0, L_inf_abs_r = 0.0;
	double L_1_rel_r = 0.0, L_2_rel_r = 0.0, L_inf_rel_r = 0.0;
	double L_1_abs   = 0.0, L_2_abs   = 0.0, L_inf_abs   = 0.0;
	double L_1_rel   = 0.0, L_2_rel   = 0.0, L_inf_rel   = 0.0;
	double* points;
	double* denoms;
	int K = 5, N = 4, M, L = N + 7, Z = 300;
	
	M = K * (N - 1) + 1;
	h = (b - a) / (M - 1);
	h_total = pow(h, N - 1);
	len = (b - a) / K;	
	delta = (b - a) / (Z - 1);
	
	MatrixXd A = MatrixXd::Zero(M, M);
	VectorXd B = VectorXd::Zero(M);
	points = new double[L * K];
	denoms = new double[N];
	
	ofstream fout("1.txt");
	fout << a << " " << b << " " << Z << " " << K * L << endl;
	
	uniform_real_distribution<double> distribution{a, a + len};
	mt19937 generator{random_device().operator ()()};
	
	cout.setf(ios::scientific);
	fout.setf(ios::scientific);

	for (int k = 0; k < K; k++){
		for (int i = 0; i < L; i++){
			points[i + k * L] = len * k + distribution(generator);
			fout << points[i + k * L] << endl;
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
	
	double tmp1 = 0.0, tmp2 = 0.0;																			// Заполнение матрицы A и вектора B
	for (int k = 0; k < K; k++){
		for (int i = k * (N - 1); i <= (k + 1) * (N - 1); i++){
			for (int j = i; j <= (k + 1) * (N - 1); j++){
				tmp1 = 0.0;
				for (int n = k * L; n < (k + 1) * L; n++){
					tmp1 += Phi_i(denoms, N, points[n], a, h, j, k) * Phi_i(denoms, N, points[n], a, h, i, k);
				}
				A(i, j) += tmp1;
				if (i != j)
					A(j, i) += tmp1;
			}
			for (int m = k * L; m < (k + 1) * L; m++){
				tmp2 += Phi_i(denoms, N, points[m], a, h, i, k) * f(points[m]);
			}
			B(i) += tmp2;			
			tmp2 = 0.0;
		}
	}
	//cout << endl << A << endl << endl << B << endl;
	VectorXd C = A.colPivHouseholderQr().solve(B);																// Решение системы
	//cout << endl << C << endl;
	
	double tmp3 = 0.0, tmp_f = 0.0;
	for (int k = 0; k < K; k++){																				// Вычисление ЭНП
		for (int i = k * Z / K; i < (k + 1) * Z / K; i++){
			for (int j = k * (N - 1); j <= (k + 1) * (N - 1); j++){
				tmp3 += C(j) * Phi_i(denoms, N, a + i * delta, a, h, j, k);
			}
			fout << tmp3 << endl;
			tmp3 = 0.0;
		}
	}
	
	////////////////////////////////////////////////////////////////////
	
	for (int k = 0; k < K; k++){
		for (int l = 0; l < L; l++){																			// Вычисление норм в "случайных" точках
			for (int j = k * (N - 1); j <= (k + 1) * (N - 1); j++){
				tmp3 += C(j) * Phi_i(denoms, N, points[l + k * L], a, h, j, k);
			}
			tmp_f = f(points[l + k * L]);
			
			L_1_abs_r   += abs(tmp_f - tmp3);
			L_2_abs_r   += (tmp_f - tmp3) * (tmp_f - tmp3);
			L_inf_abs_r  = max(L_inf_abs_r, abs(tmp_f - tmp3));
			
			L_1_f   += abs(tmp_f);
			L_2_f   += tmp_f * tmp_f;
			L_inf_f  = max(L_inf_f, abs(tmp_f));
			
			tmp3 = 0.0;
		}
	}
	
	L_2_abs_r = sqrt(L_2_abs_r);
	L_2_f     = sqrt(L_2_f);
	
	L_1_rel_r   = L_1_abs_r / L_1_f;
	L_2_rel_r   = L_2_abs_r / L_2_f;
	L_inf_rel_r = L_inf_abs_r / L_inf_f;
	
	L_1_f   = 0.0;
	L_2_f   = 0.0;
	L_inf_f = 0.0;
	
	cout << "Random points |      L1      |      L2      |     Linf     |" << endl;
	cout << "Absolute      |" << setw(14) << L_1_abs_r << "|" << setw(14)  << L_2_abs_r << "|" << setw(14) << L_inf_abs_r << "|" << endl;
	cout << "Relative      |" << setw(14) << L_1_rel_r << "|" << setw(14)  << L_2_rel_r << "|" << setw(14) << L_inf_rel_r << "|" << endl;
	cout << endl;
	
	//////////////////////////////////////////////////////////////////
	
	for (int k = 0; k < K; k++){																		// Вычисление норм на сетке с шагом h/100
		for (int i = k * (100 * (M - 1) + 1) / K; i < (k + 1) * (100 * (M - 1) + 1) / K; i++){
			for (int j = k * (N - 1); j <= (k + 1) * (N - 1); j++){
				tmp3 += C(j) * Phi_i(denoms, N, a + i * h / 100, a, h, j, k);
			}
			tmp_f = f(a + i * h / 100);
			
			L_1_abs   += abs(tmp_f - tmp3);
			L_2_abs   += (tmp_f - tmp3) * (tmp_f - tmp3);
			L_inf_abs  = max(L_inf_abs, abs(tmp_f - tmp3));
			
			L_1_f   += abs(tmp_f);
			L_2_f   += tmp_f * tmp_f;
			L_inf_f  = max(L_inf_f, abs(tmp_f));
			
			tmp3 = 0.0;
		}
	}
	
	L_2_abs = sqrt(L_2_abs);
	L_2_f   = sqrt(L_2_f);
	
	L_1_rel   = L_1_abs / L_1_f;
	L_2_rel   = L_2_abs / L_2_f;
	L_inf_rel = L_inf_abs / L_inf_f;
	
	cout << "h/100 points  |      L1      |      L2      |     Linf     |" << endl;
	cout << "Absolute      |" << setw(14) << L_1_abs << "|" << setw(14) << L_2_abs << "|" << setw(14) << L_inf_abs << "|" << endl;
	cout << "Relative      |" << setw(14) << L_1_rel << "|" << setw(14) << L_2_rel << "|" << setw(14) << L_inf_rel << "|" << endl;
	
	
	fout.close();
	
	delete [] points;
	delete [] denoms;
	
	system("py task3.py");
	
	return 0;
}
