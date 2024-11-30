#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;


double f(double x, double t){
	return t * cos(3.14 * x);
}


double Psi_0_f(double t){
	return t * t;
}


double Psi_1_f(double t){
	return 1. - t;
}


double Phi_f(double x){
	return x;
}


double** Explicit(double h, double tau, int M, int N){							// Явная схема
	double **U = new double*[N];
	
	for (int i = 0; i < N; i++){
		U[i] = new double[M];
	}
	
	for (int j = 0; j < M; j++){
		U[0][j] = Phi_f(j * h);
	}
	
	for (int i = 0; i < N; i++){
		U[i][0]     = Psi_0_f(i * tau);
		U[i][M - 1] = Psi_1_f(i * tau);
	}
	
	for (int i = 0; i < N - 1; i++){
		for (int j = 1; j < M - 1; j++){
			U[i + 1][j] = U[i][j] + tau/pow(h, 2) * (U[i][j - 1] - 2. * U[i][j] + U[i][j + 1]) + f(h * j, tau * i) * tau;
		}
	}
	
	return U;
}


double** Implicit(double h, double tau, int M, int N, double sigma){			// Неявная схема с весами
	double alpha, beta, gamma, denom;
	double **U = new double*[N];
	double  *A = new double[M - 3];
	double  *B = new double[M - 3];
	
	for (int i = 0; i < N; i++){
		U[i] = new double[M];
	}
	
	for (int i = 0; i < M; i++){
		U[0][i] = Phi_f(i * h);
	}
	
	for (int i = 0; i < N; i++){
		U[i][0]     = Psi_0_f(i * tau);
		U[i][M - 1] = Psi_1_f(i * tau);
	}
	
	alpha = 1./tau + (2. * sigma)/pow(h, 2);
	beta  = - sigma/pow(h, 2);
	gamma = - sigma/pow(h, 2);
	
	A[0] = - beta/alpha;
	for (int j = 0; j < M - 4; j++){
		denom = gamma * A[j] + alpha;
			
		A[j + 1] = - beta/denom;
	}
	
	for (int i = 0; i < N - 1; i++){
		B[0] = ((1. - sigma) * (U[i][2] - 2. * U[i][1] + U[i][0])/pow(h, 2) + f(h, tau * i) * (1. - sigma) + 
				f(h, tau * (i + 1)) * sigma + sigma/pow(h, 2) * U[i + 1][0] + U[i][1]/tau)/alpha;
	
		for (int j = 0; j < M - 4; j++){
			denom = gamma * A[j] + alpha;
			
			B[j + 1] = ((1. - sigma) * (U[i][j + 3] - 2. * U[i][j + 2] + U[i][j + 1])/pow(h, 2) + f(h * (j + 2), tau * i) * (1. - sigma) + 
						f(h * (j + 2), tau * (i + 1)) * sigma + U[i][j + 2]/tau - gamma * B[j])/denom;
		}
		
		U[i + 1][M - 2] = - (B[M - 4] - ((1. - sigma) * (U[i][M - 1] - 2. * U[i][M - 2] + U[i][M - 3])/pow(h, 2) + f(h * (M - 2), tau * i) * (1. - sigma) + 
							f(h * (M - 2), tau * (i + 1)) * sigma + sigma/pow(h, 2) * U[i + 1][M - 1] + U[i][M - 2]/tau)/gamma) / (A[M - 4] + alpha/gamma);
		for (int j = M - 3; j > 0; j--){
			U[i + 1][j] = A[j - 1] * U[i + 1][j + 1] + B[j - 1];
		}
	}
	
	delete [] A;
	delete [] B;
	
	return U;
}


void error(double** Y_0, double** Y_1, int N, int M, int k){
	double L_1_y_abs = 0; 		// L_1 абсолютное отклонение h относительно h/2
	double L_2_y_abs = 0; 		// L_2 абсолютное отклонение h относительно h/2
	double L_i_y_abs = 0;		// L_i абсолютное отклонение h относительно h/2
	double L_1_y = 0;			// L_1 норма решения y на сетке с шагом h
	double L_2_y = 0;			// L_2 норма решения y на сетке с шагом h
	double L_i_y = 0;			// L_i норма решения y на сетке с шагом h
	
	//ofstream fout("errors.txt");
	//fout.setf(ios::scientific);
	
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			L_1_y 	  += abs(Y_0[i][j]);
			L_1_y_abs += abs(Y_0[i][j] - Y_1[k * i][2 * j]);
			
			L_2_y 	  += Y_0[i][j] * Y_0[i][j];
			L_2_y_abs += (Y_0[i][j] - Y_1[k * i][2 * j]) * (Y_0[i][j] - Y_1[k * i][2 * j]);
			
			L_i_y 	  = max(L_i_y, Y_0[i][j]);
			L_i_y_abs = max(L_i_y_abs, abs(Y_0[i][j] - Y_1[k * i][2 * j]));
		}
	}
	
	L_2_y_abs = sqrt(L_2_y_abs);
	L_2_y 	  = sqrt(L_2_y);
	
	
	cout.setf(ios::scientific);
	
	cout << "              |-----------------------------|" << endl;
	cout << "              |   Absolute   |   Relative   |" << endl;
	cout << "              |  h/2, tau/" << k << "  |  h/2, tau/" << k << "  |" << endl;
	cout << "L_1 error     |" << setw(14) << L_1_y_abs << '|' << setw(14) << L_1_y_abs / L_1_y << '|' << endl;
	cout << "L_2 error     |" << setw(14) << L_2_y_abs << '|' << setw(14) << L_2_y_abs / L_2_y << '|' << endl;
	cout << "L_i error     |" << setw(14) << L_i_y_abs << '|' << setw(14) << L_i_y_abs / L_i_y << '|' << endl;
	cout << "              |-----------------------------|" << endl;
	
}


int main(void){
	int M_1 = 15, N_1, M_2 = 15, N_2 = 55;
	double h_1, tau_1, h_2, tau_2, sigma = 0.5;
	double **U_1;
	double **U_2;
	
	ofstream fout1("1.txt");
	ofstream fout2("2.txt");
	fout1.setf(ios::scientific);
	fout2.setf(ios::scientific);
	
	N_1 = 2 * (M_1 - 1) * (M_1 - 1) + 1;
	//N_2 = 2 * (M_2 - 1) * (M_2 - 1) + 1;
	
	h_1   = 1./(M_1 - 1);
	tau_1 = 1./(N_1 - 1);
	h_2   = 1./(M_2 - 1);
	tau_2 = 1./(N_2 - 1);
	
	
	U_1 = Explicit(h_1,   tau_1,   M_1,         N_1);
	U_2 = Explicit(h_1/2, tau_1/4, 2 * M_1 - 1, 4 * N_1 - 3);
	for (int i = 0; i < N_1; i++){
		for (int j = 0; j < M_1; j++){
			fout1 << setw(14) << U_1[i][j] << " ";
		}
		fout1 << endl;
	}
	cout << "              |           Explicit          |" << endl;
	error(U_1, U_2, N_1, M_1, 4);
	
	for (int i = 0; i < N_1; i++){
		delete [] U_1[i];
		delete [] U_2[i];
	}
	delete [] U_1;
	delete [] U_2;
	
	
	U_1 = Implicit(h_2,   tau_2,   M_2,         N_2,         sigma);
	U_2 = Implicit(h_2/2, tau_2/2, 2 * M_2 - 1, 2 * N_2 - 1, sigma);
	for (int i = 0; i < N_2; i++){
		for (int j = 0; j < M_2; j++){
			fout2 << setw(14) << U_1[i][j] << " ";
		}
		fout2 << endl;
	}
	cout << "              |           Implicit          |" << endl;
	error(U_1, U_2, N_2, M_2, 2);
	
	
	for (int i = 0; i < N_2; i++){
		delete [] U_1[i];
		delete [] U_2[i];
	}
	delete [] U_1;
	delete [] U_2;
	
	fout1.close();
	fout2.close();
	
	system("py task5.py");
	
	return 0;
}