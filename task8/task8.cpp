#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;



double f(double x, double y, double z){
	return z * sqrt(x) + y * pow(2.718, x + 1) + cos(x);
}


double** Euler(double h, int N, double a, double q_0, double q_1){
	double x, y, z;
	double x_p = a;
	double y_p = q_0;
	double z_p = q_1;
	double **Y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		Y[i] = new double[N];
	}
	
	Y[0][0] = y_p;
	Y[1][0] = z_p;
	
	for (int i = 1; i < N; i++){												// Метод Эйлера
		x = x_p + h;
		y = y_p + h * z_p;
		z = z_p + h * f(x_p, y_p, z_p);		
		
		Y[0][i] = y;
		Y[1][i] = z;
		
		x_p = x;
		y_p = y;
		z_p = z;
	}
		
	return Y;
}


double** Euler_with_recount(double h, int N, double a, double q_0, double q_1){
	double x, y, z, y_predict, z_predict;
	double x_p = a;
	double y_p = q_0;
	double z_p = q_1;
	double **Y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		Y[i] = new double[N];
	}
	
	Y[0][0] = y_p;
	Y[1][0] = z_p;
	
	for (int i = 1; i < N; i++){												// Метод Эйлера с пересчетом
		y_predict = y_p + h * z_p;
		z_predict = z_p + h * f(x_p, y_p, z_p);
		
		x = x_p + h;
		y = y_p + h * z_predict;
		z = z_p + h * f(x_p, y_predict, z_predict);
		
		Y[0][i] = y;
		Y[1][i] = z;
		
		x_p = x;
		y_p = y;
		z_p = z;
	}
	
	return Y;
}


double** Runge_Kutta_2_order(double h, int N, double a, double q_0, double q_1){
	double k_1, k_2, k_11, k_21, x, y, z;
	double x_p = a;
	double y_p = q_0;
	double z_p = q_1;
	double **Y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		Y[i] = new double[N];
	}
	
	Y[0][0] = y_p;
	Y[1][0] = z_p;
	
	for (int i = 1; i < N; i++){												// Метод Рунге-Кутта 2 порядка
		k_1  = f(x_p,     y_p,           z_p);
		k_2  = f(x_p + h, y_p + k_1 * h, z_p + k_1 * h);
		k_11 = z_p;
		k_21 = z_p + h * k_11;
		
		x = x_p + h;
		y = y_p + (k_11 + k_21) * h/2;
		z = z_p + (k_1  + k_2)  * h/2;
		
		Y[0][i] = y;
		Y[1][i] = z;
		
		x_p = x;
		y_p = y;
		z_p = z;
	}
	
	return Y;
}


double** Runge_Kutta_4_order(double h, int N, double a, double q_0, double q_1){
	double k_1, k_2, k_3, k_4, k_11, k_21, k_31, k_41, x, y, z;
	double x_p = a;
	double y_p = q_0;
	double z_p = q_1;
	double **Y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		Y[i] = new double[N];
	}
	
	Y[0][0] = y_p;
	Y[1][0] = z_p;
	
	for (int i = 1; i < N; i++){												// Метод Рунге-Кутта 4 порядка
		k_1 = f(x_p, y_p, z_p);
		k_2 = f(x_p + h/2, y_p + k_1 * h/2, z_p + k_1 * h/2);
		k_3 = f(x_p + h/2, y_p + k_2 * h/2, z_p + k_2 * h/2);
		k_4 = f(x_p + h,   y_p + k_3 * h,   z_p + k_3 * h);
		
		k_11 = z_p;
		k_21 = z_p + h * k_11/2;
		k_31 = z_p + h * k_21/2;
		k_41 = z_p + h * k_31;
		
		x = x_p + h;
		y = y_p + (k_11 + 2 * k_21 + 2 * k_31 + k_41) * h/6;
		z = z_p + (k_1  + 2 * k_2  + 2 * k_3  + k_4)  * h/6;
		
		Y[0][i] = y;
		Y[1][i] = z;
		
		x_p = x;
		y_p = y;
		z_p = z;
	}
	
	return Y;
}


double** Adams(double h, int N, double a, double q_0, double q_1){
	double K_1, K_2, K_3, K_11, K_21, K_31, k_1, k_2, k_3, k_4, k_11, k_21, k_31, k_41, x, y, z;
	double x_n, y_n, z_n, x_n_1, y_n_1, z_n_1, x_n_2, y_n_2, z_n_2;
	double **Y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		Y[i] = new double[N];
	}
	
	x_n_2 = a;
	y_n_2 = q_0;
	z_n_2 = q_1;
	
	Y[0][0] = y_n_2;
	Y[1][0] = z_n_2;
	
	k_1 = f(x_n_2, y_n_2, z_n_2);
	k_2 = f(x_n_2 + h/2, y_n_2 + k_1 * h/2, z_n_2 + k_1 * h/2);
	k_3 = f(x_n_2 + h/2, y_n_2 + k_2 * h/2, z_n_2 + k_2 * h/2);
	k_4 = f(x_n_2 + h,   y_n_2 + k_3 * h,   z_n_2 + k_3 * h);
	
	k_11 = z_n_2;
	k_21 = z_n_2 + h * k_11/2;
	k_31 = z_n_2 + h * k_21/2;
	k_41 = z_n_2 + h * k_31;
	
	x_n_1 = x_n_2 + h;
	y_n_1 = y_n_2 + (k_11 + 2 * k_21 + 2 * k_31 + k_41) * h/6;
	z_n_1 = z_n_2 + (k_1  + 2 * k_2  + 2 * k_3  + k_4)  * h/6;
	
	Y[0][1] = y_n_1;
	Y[1][1] = z_n_1;
	
	k_1 = f(x_n_1, y_n_1, z_n_1);
	k_2 = f(x_n_1 + h/2, y_n_1 + k_1 * h/2, z_n_1 + k_1 * h/2);
	k_3 = f(x_n_1 + h/2, y_n_1 + k_2 * h/2, z_n_1 + k_2 * h/2);
	k_4 = f(x_n_1 + h,   y_n_1 + k_3 * h,   z_n_1 + k_3 * h);
	
	k_11 = z_n_1;
	k_21 = z_n_1 + h * k_11/2;
	k_31 = z_n_1 + h * k_21/2;
	k_41 = z_n_1 + h * k_31;
	
	x_n = x_n_1 + h;
	y_n = y_n_1 + (k_11 + 2 * k_21 + 2 * k_31 + k_41) * h/6;
	z_n = z_n_1 + (k_1  + 2 * k_2  + 2 * k_3  + k_4)  * h/6;
	
	Y[0][2] = y_n;
	Y[1][2] = z_n;
	
	for (int i = 3; i < N; i++){
		K_1 = f(x_n,   y_n,   z_n);
		K_2 = f(x_n_1, y_n_1, z_n_1);
		K_3 = f(x_n_2, y_n_2, z_n_2);
		
		K_11 = z_n;
		K_21 = z_n_1;
		K_31 = z_n_2;		
		
		x = x_n + h;
		y = y_n + h * (23./12. * K_11 - 16./12. * K_21 + 5./12. * K_31);
		z = z_n + h * (23./12. * K_1  - 16./12. * K_2  + 5./12. * K_3);
		
		Y[0][i] = y;
		Y[1][i] = z;
		
		x_n_2 = x_n_1;
		y_n_2 = y_n_1;
		z_n_2 = z_n_1;
		
		x_n_1 = x_n;
		y_n_1 = y_n;
		z_n_1 = z_n;
		
		x_n = x;
		y_n = y;
		z_n = z;
	}
	
	return Y;
}


void error(double** Y_0, double** Y_1, double** Y_2, int N){
	double L_1_1_y_abs = 0, L_1_1_z_abs = 0; 	// L_1 абсолютное отклонение y и y' h относительно h/2
	double L_1_2_y_abs = 0, L_1_2_z_abs = 0; 	// L_1 абсолютное отклонение y и y' h относительно 2h
	double L_2_1_y_abs = 0, L_2_1_z_abs = 0; 	// L_2 абсолютное отклонение y и y' h относительно h/2
	double L_2_2_y_abs = 0, L_2_2_z_abs = 0; 	// L_2 абсолютное отклонение y и y' h относительно 2h
	double L_i_1_y_abs = 0, L_i_1_z_abs = 0;	// L_i абсолютное отклонение y и y' h относительно h/2
	double L_i_2_y_abs = 0, L_i_2_z_abs = 0;	// L_i абсолютное отклонение y и y' h относительно 2h
	double L_1_y = 0, L_1_z = 0;				// L_1 норма решения y и y' на сетке с шагом h
	double L_2_y = 0, L_2_z = 0;				// L_2 норма решения y и y' на сетке с шагом h
	double L_i_y = 0, L_i_z = 0;				// L_i норма решения y и y'на сетке с шагом h
	
	//ofstream fout("errors.txt");
	//fout.setf(ios::scientific);
	
	
	for (int i = 0; i < N; i++){
		L_1_y 		+= abs(Y_0[0][i]);
		L_1_z 		+= abs(Y_0[1][i]);
		L_1_1_y_abs += abs(Y_0[0][i] - Y_1[0][2 * i]);
		L_1_1_z_abs += abs(Y_0[1][i] - Y_1[1][2 * i]);
		if (i % 2 == 0){
			L_1_2_y_abs += abs(Y_0[0][i] - Y_2[0][i/2]);
			L_1_2_z_abs += abs(Y_0[1][i] - Y_2[1][i/2]);
		}
		
		L_2_y 		+= Y_0[0][i] * Y_0[0][i];
		L_2_z 		+= Y_0[1][i] * Y_0[1][i];
		L_2_1_y_abs += (Y_0[0][i] - Y_1[0][2 * i]) * (Y_0[0][i] - Y_1[0][2 * i]);
		L_2_1_z_abs += (Y_0[1][i] - Y_1[1][2 * i]) * (Y_0[1][i] - Y_1[1][2 * i]);
		if (i % 2 == 0){
			L_2_2_y_abs += (Y_0[0][i] - Y_2[0][i/2]) * (Y_0[0][i] - Y_2[0][i/2]);
			L_2_2_z_abs += (Y_0[1][i] - Y_2[1][i/2]) * (Y_0[1][i] - Y_2[1][i/2]);
		}
		
		L_i_y 		= max(L_i_y, Y_0[0][i]);
		L_i_z 		= max(L_i_z, Y_0[1][i]);
		L_i_1_y_abs = max(L_i_1_y_abs, abs(Y_0[0][i] - Y_1[0][2 * i]));
		L_i_1_z_abs = max(L_i_1_z_abs, abs(Y_0[1][i] - Y_1[1][2 * i]));
		if (i % 2 == 0){
			L_i_2_y_abs = max(L_i_2_y_abs, abs(Y_0[0][i] - Y_2[0][i/2]));
			L_i_2_z_abs = max(L_i_2_z_abs, abs(Y_0[1][i] - Y_2[1][i/2]));
		}
		//cout << L_2_1_y_abs << " " << L_2_1_z_abs << endl;
	}
	
	L_2_1_y_abs = sqrt(L_2_1_y_abs);
	L_2_2_y_abs = sqrt(L_2_2_y_abs);
	L_2_y 	    = sqrt(L_2_y);
	L_2_1_z_abs = sqrt(L_2_1_z_abs);
	L_2_2_z_abs = sqrt(L_2_2_z_abs);
	L_2_z	    = sqrt(L_2_z);
	
	//cout << L_i_y << " " << L_i_z;
	
	cout.setf(ios::scientific);
	
	cout << "              |-----------------------------------------------------------|" << endl;
	cout << "              |                             y                             |" << endl;
	cout << "              |          Absolute           |          Relative           |" << endl;
	cout << "              |      h/2     |      2*h     |      h/2     |      2*h     |" << endl;
	cout << "L_1 error     |" << setw(14) << L_1_1_y_abs << '|' << setw(14) << L_1_2_y_abs << '|';
	cout << setw(14) << L_1_1_y_abs / L_1_y << '|' << setw(14) << L_1_2_y_abs / L_1_y << '|' << endl;
	cout << "L_2 error     |" << setw(14) << L_2_1_y_abs << '|' << setw(14) << L_2_2_y_abs << '|';
	cout << setw(14) << L_2_1_y_abs / L_2_y << '|' << setw(14) << L_2_2_y_abs / L_2_y << '|' << endl;
	cout << "L_i error     |" << setw(14) << L_i_1_y_abs << '|' << setw(14) << L_i_2_y_abs << '|';
	cout << setw(14) << L_i_1_y_abs / L_i_y << '|' << setw(14) << L_i_2_y_abs / L_i_y << '|' << endl;
	cout << "              |-----------------------------------------------------------|" << endl;
	cout << "              |                             y'                            |" << endl;
	cout << "              |          Absolute           |          Relative           |" << endl;
	cout << "              |      h/2     |      2*h     |      h/2     |      2*h     |" << endl;
	cout << "L_1 error     |" << setw(14) << L_1_1_z_abs << '|' << setw(14) << L_1_2_z_abs << '|';
	cout << setw(14) << L_1_1_z_abs / L_1_z << '|' << setw(14) << L_1_2_z_abs / L_1_z << '|' << endl;
	cout << "L_2 error     |" << setw(14) << L_2_1_z_abs << '|' << setw(14) << L_2_2_z_abs << '|';
	cout << setw(14) << L_2_1_z_abs / L_2_z << '|' << setw(14) << L_2_2_z_abs / L_2_z << '|' << endl;
	cout << "L_i error     |" << setw(14) << L_i_1_z_abs << '|' << setw(14) << L_i_2_z_abs << '|';
	cout << setw(14) << L_i_1_z_abs / L_i_z << '|' << setw(14) << L_i_2_z_abs / L_i_z << '|' << endl;
	cout << "              |-----------------------------------------------------------|" << endl;
	
}


int main(void){
	double a = 0., b = 1., q_0 = 0., q_1 = 1., h_1, h_2, h_3, h_4, h_5;
	int N_1 = 55, N_2 = 55, N_3 = 55, N_4 = 55, N_5 = 55;
	double **Y_0;
	double **Y_1;
	double **Y_2;
	
	
	h_1 = (b - a) / (N_1 - 1);
	h_2 = (b - a) / (N_2 - 1);
	h_3 = (b - a) / (N_3 - 1);
	h_4 = (b - a) / (N_4 - 1);
	h_5 = (b - a) / (N_5 - 1);
	
	ofstream fout1("1.txt");
	ofstream fout2("2.txt");
	ofstream fout3("3.txt");
	ofstream fout4("4.txt");
	ofstream fout5("5.txt");
	fout1.setf(ios::scientific);
	fout2.setf(ios::scientific);
	fout3.setf(ios::scientific);
	fout4.setf(ios::scientific);
	fout5.setf(ios::scientific);
	
	
	Y_0 = Euler(h_1,     N_1,             a, q_0, q_1);
	Y_1 = Euler(h_1/2,   2 * N_1 - 1,     a, q_0, q_1);
	Y_2 = Euler(2 * h_1, (N_1 - 1)/2 + 1, a, q_0, q_1);
	
	for (int i = 0; i < N_1; i++){
		fout1 << a + h_1 * i << " " << Y_0[0][i] << " " << Y_0[1][i] << endl;
	}
	cout << "              |                           Euler                           |" << endl;
	error(Y_0, Y_1, Y_2, N_1);
	
	for (int i = 0; i < 2; i++){
		delete [] Y_0[i];
		delete [] Y_1[i];
		delete [] Y_2[i];
	}
	delete [] Y_0;
	delete [] Y_1;
	delete [] Y_2;
	
	
	Y_0 = Euler_with_recount(h_2,     N_2,             a, q_0, q_1);
	Y_1 = Euler_with_recount(h_2/2,   2 * N_2 - 1,     a, q_0, q_1);
	Y_2 = Euler_with_recount(2 * h_2, (N_2 - 1)/2 + 1, a, q_0, q_1);
	
	for (int i = 0; i < N_2; i++){
		fout2 << a + h_2 * i << " " << Y_0[0][i] << " " << Y_0[1][i] << endl;
	}
	cout << "              |                    Euler with recount                     |" << endl;
	error(Y_0, Y_1, Y_2, N_2);
	
	for (int i = 0; i < 2; i++){
		delete [] Y_0[i];
		delete [] Y_1[i];
		delete [] Y_2[i];
	}
	delete [] Y_0;
	delete [] Y_1;
	delete [] Y_2;
	
	
	Y_0 = Runge_Kutta_2_order(h_3,     N_3,             a, q_0, q_1);
	Y_1 = Runge_Kutta_2_order(h_3/2,   2 * N_3 - 1,     a, q_0, q_1);
	Y_2 = Runge_Kutta_2_order(2 * h_3, (N_3 - 1)/2 + 1, a, q_0, q_1);
	
	for (int i = 0; i < N_3; i++){
		fout3 << a + h_3 * i << " " << Y_0[0][i] << " " << Y_0[1][i] << endl;
	}
	cout << "              |                   Runge-Kutta 2nd order                   |" << endl;
	error(Y_0, Y_1, Y_2, N_3);
	
	for (int i = 0; i < 2; i++){
		delete [] Y_0[i];
		delete [] Y_1[i];
		delete [] Y_2[i];
	}
	delete [] Y_0;
	delete [] Y_1;
	delete [] Y_2;
	
	
	Y_0 = Runge_Kutta_4_order(h_4,     N_4,             a, q_0, q_1);
	Y_1 = Runge_Kutta_4_order(h_4/2,   2 * N_4 - 1,     a, q_0, q_1);
	Y_2 = Runge_Kutta_4_order(2 * h_4, (N_4 - 1)/2 + 1, a, q_0, q_1);
	
	for (int i = 0; i < N_4; i++){
		fout4 << a + h_4 * i << " " << Y_0[0][i] << " " << Y_0[1][i] << endl;
	}
	cout << "              |                   Runge-Kutta 4th order                   |" << endl;
	error(Y_0, Y_1, Y_2, N_4);
	
	for (int i = 0; i < 2; i++){
		delete [] Y_0[i];
		delete [] Y_1[i];
		delete [] Y_2[i];
	}
	delete [] Y_0;
	delete [] Y_1;
	delete [] Y_2;
	
	
	Y_0 = Adams(h_5,     N_5,             a, q_0, q_1);
	Y_1 = Adams(h_5/2,   2 * N_5 - 1,     a, q_0, q_1);
	Y_2 = Adams(2 * h_5, (N_5 - 1)/2 + 1, a, q_0, q_1);
	
	for (int i = 0; i < N_5; i++){
		fout5 << a + h_5 * i << " " << Y_0[0][i] << " " << Y_0[1][i] << endl;
	}
	cout << "              |                           Adams                           |" << endl;
	error(Y_0, Y_1, Y_2, N_5);
	
	for (int i = 0; i < 2; i++){
		delete [] Y_0[i];
		delete [] Y_1[i];
		delete [] Y_2[i];
	}
	delete [] Y_0;
	delete [] Y_1;
	delete [] Y_2;
	
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
	
	system("py task3.py");
	
	return 0;
}

