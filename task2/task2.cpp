#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;


double f(double x){
	//return 1.0 / (1.0 + x * x);
	return x * sin(x);
	//return 0.15 *  x * x * x * x + 2 * x * x + 1;
	//return 2 * x;
}



double Lagrange(double* denoms, int size, double x, double length, double step, double a, int section){				// Многочлен Лагранжа
	double value_x = 0.0, num = 1.0;
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			if (j != i){
				num *= (x - a - length * section - step * j);
			}
		}
		value_x += f(a + length * section + step * i) * num / (denoms[i]);
		num = 1.0;
	}
	return value_x;
}


int main(void){
	double a = -3.0, b = 3.0, den = 1.0, h, L, delta, L_2 = 0.0, L_1 = 0.0, L_inf = 0, tmp1 = 0, tmp2 = 0, h_total;	
	double L_1_rel = 0.0, L_2_rel = 0.0, L_inf_rel = 0.0, L_1_f = 0.0, L_2_f = 0.0, L_inf_f = 0.0;
	double *denoms;
	int N = 4, K = 10, Z = 300, M;
	
	M = K * (N - 1) + 1;
	h = (b - a) / (M - 1);
	h_total = pow(h, N - 1);
	L = (b - a) / K;
	delta = (b - a) / (Z - 1);
	denoms = new double[N];
	
	ofstream fout("1.txt");
	fout << a << " " << b << " " << Z << " " << M << endl;
	
	for (int j = 0; j < N; j++){																					// Заполнение массива знаменателей
		for (int k = 0; k < N; k++){
			if (k != j)
				den *= (j - k);
		}
		denoms[j] = den * h_total;
		den = 1.0;
	}
	
	for (int k = 0; k < K; k++){																				// Вычисление полинома
		for (int i = k * Z / K; i < (k + 1) * Z / K; i++){
			fout << Lagrange(denoms, N, a + i * delta, L, h, a, k) << endl;
		}
	}
	
	for (int k = 0; k < K; k++){																		// Вычисление норм
		for (int i = k * (100 * (M - 1) + 1) / K; i < (k + 1) * (100 * (M - 1) + 1) / K; i++){
			tmp1 = Lagrange(denoms, N, a + i * h / 100, L, h, a, k);
			tmp2 = f(a + i * h / 100);
			L_2   += (tmp1 - tmp2) * (tmp1 - tmp2);
			L_2_f += tmp2 * tmp2;
			L_1_f += abs(tmp2);
			L_1   += abs((tmp1 - tmp2));
			if (abs(tmp1 - tmp2) > L_inf)
				L_inf = abs(tmp1 - tmp2);
			if (abs(tmp2) > L_inf_f)
				L_inf_f = abs(tmp2);
		}
	}
	L_2       = sqrt(L_2);
	L_2_f     = sqrt(L_2_f);
	L_2_rel   = L_2 / L_2_f;
	L_1_rel   = L_1 / L_1_f;
	L_inf_rel = L_inf / L_inf_f;
	
	
	cout.setf(ios::scientific);
	//cout << "L1-norm(absolute):" << L_1     << endl << "L2-norm(absolute):" << L_2     << endl << "Linf-norm(absolute):" << L_inf << endl;
	//cout << "L1-norm(relative):" << L_1_rel << endl << "L2-norm(relative):" << L_2_rel << endl << "Linf-norm(relative):" << L_inf_rel;
	
	cout << "          |      L1      |      L2      |     Linf     |" << endl;
	cout << "Absolute  |" << setw(14) << L_1     << "|" << setw(14) << L_2     << "|" << setw(14) << L_inf     << "|" << endl;
	cout << "Relative  |" << setw(14) << L_1_rel << "|" << setw(14) << L_2_rel << "|" << setw(14) << L_inf_rel << "|" << endl;
	
	fout.close();
	
	delete [] denoms;
	system("py task2.py");
	
	return 0;
}