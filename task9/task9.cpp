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


double f_h(double x, double y, double z){
	return z * sqrt(x) + y * pow(2.718, x + 1);
}


double f_nl(double x, double y, double z){
	return pow(z, 2) * sqrt(x) + pow(y, 2) * pow(2.718, x + 1) + cos(x);
}


double G(double x, double y, double u, double Y, double U){
	return 2 * y * pow(2.718, x + 1) * Y + 2 * u * sqrt(x) * U;
}


double RK_4_YU(double h, int N, double a, double q_0, double q_1, double *y, double *z){
	double k_1, k_2, k_3, k_4, k_11, k_21, k_31, k_41, x, Y, U;
	double x_p = a;
	double Y_p = q_0;
	double U_p = q_1;
	
	for (int i = 1; i < N; i++){												// Метод Рунге-Кутта 4 порядка
		k_1 = G(x_p,       y[2 * i - 2], z[2 * i - 2], Y_p,             U_p);
		k_2 = G(x_p + h/2, y[2 * i - 1], z[2 * i - 1], Y_p + k_1 * h/2, U_p + k_1 * h/2);
		k_3 = G(x_p + h/2, y[2 * i - 1], z[2 * i - 1], Y_p + k_2 * h/2, U_p + k_2 * h/2);
		k_4 = G(x_p + h,   y[2 * i],     z[2 * i],     Y_p + k_3 * h,   U_p + k_3 * h);
		
		k_11 = U_p;
		k_21 = U_p + h * k_11/2;
		k_31 = U_p + h * k_21/2;
		k_41 = U_p + h * k_31;
		
		x = x_p + h;
		Y = Y_p + (k_11 + 2 * k_21 + 2 * k_31 + k_41) * h/6;
		U = U_p + (k_1  + 2 * k_2  + 2 * k_3  + k_4) * h/6;
		
		x_p = x;
		Y_p = Y;
		U_p = U;
	}
	
	return Y;
}


double** Runge_Kutta_4_order_full(double (*F)(double, double, double), double h, int N, double a, double q_0, double q_1){
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
		k_1 = F(x_p, y_p, z_p);
		k_2 = F(x_p + h/2, y_p + k_1 * h/2, z_p + k_1 * h/2);
		k_3 = F(x_p + h/2, y_p + k_2 * h/2, z_p + k_2 * h/2);
		k_4 = F(x_p + h,   y_p + k_3 * h,   z_p + k_3 * h);
		
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


double* Runge_Kutta_4_order(double (*F)(double, double, double), double h, int N, double a, double q_0, double q_1){
	double k_1, k_2, k_3, k_4, k_11, k_21, k_31, k_41, x, y, z;
	double x_p = a;
	double y_p = q_0;
	double z_p = q_1;
	double *Y = new double[N];
	
	Y[0] = y_p;
	
	for (int i = 1; i < N; i++){												// Метод Рунге-Кутта 4 порядка
		k_1 = F(x_p, y_p, z_p);
		k_2 = F(x_p + h/2, y_p + k_1 * h/2, z_p + k_1 * h/2);
		k_3 = F(x_p + h/2, y_p + k_2 * h/2, z_p + k_2 * h/2);
		k_4 = F(x_p + h,   y_p + k_3 * h,   z_p + k_3 * h);
		
		k_11 = z_p;
		k_21 = z_p + h * k_11/2;
		k_31 = z_p + h * k_21/2;
		k_41 = z_p + h * k_31;
		
		x = x_p + h;
		y = y_p + (k_11 + 2 * k_21 + 2 * k_31 + k_41) * h/6;
		z = z_p + (k_1  + 2 * k_2  + 2 * k_3  + k_4)  * h/6;
		
		Y[i] = y;
		
		x_p = x;
		y_p = y;
		z_p = z;
	}
	
	return Y;
}


double* Shooting(double h, int N, double a, double q_a, double q_b){
	double *Y   = new double[N];
	double *Y_h = new double[N];
	double *Y_b = new double[N];
	
	Y   = Runge_Kutta_4_order(f,   h, N, a, q_a, 0.);
	Y_h = Runge_Kutta_4_order(f_h, h, N, a, 0. , 1.);
	for (int i = 0; i < N; i++){
		Y_b[i] = Y[i] + (q_b - Y[N - 1]) / Y_h[N - 1] * Y_h[i];
	}
	//cout << (q_b - Y[N - 1]) / Y_h[N - 1];
	
	delete [] Y;
	delete [] Y_h;
	
	return Y_b;
}


double* Tridiagonal(double h, int N, double a, double q_a, double q_b){
	double *Y = new double[N];
	double *A = new double[N - 1];
	double *B = new double[N - 1];
	double denom;
	
	A[0] = 0;
	B[0] = q_a;
		
	for (int i = 0; i < N - 2; i++){
		denom = (1/pow(h, 2) + sqrt(a + i * h)/(2 * h)) * A[i] + (- 2/pow(h, 2) - pow(2.718, a + i * h + 1));
		
		A[i + 1] = - (1/pow(h, 2) - sqrt(a + i * h)/(2 * h)) / denom;
		B[i + 1] = (cos(a + i * h) - (1/pow(h, 2) + sqrt(a + i * h)/(2 * h)) * B[i]) / denom;
	}
		
	Y[0]     = q_a;
	Y[N - 1] = q_b;
	
	for (int i = N - 2; i > 0; i--){
		Y[i] = A[i] * Y[i + 1] + B[i];
	}
	
	delete [] A;
	delete [] B;
		
	return Y;
}


double* Newton(double h, int N, double a, double q_a, double q_b){
	int k = 1;
	double s = 0.5, Y, eps = 1e-4, B;
	double *y_f = new double[N];
	double **y = new double*[2];
	
	for (int i = 0; i < 2; i++){
		y[i] = new double[N];
	}
	
	y = Runge_Kutta_4_order_full(f_nl, h, N, a, q_a, s);
	Y = RK_4_YU(2 * h, (N - 1)/2 + 1, a, 0, 1, y[0], y[1]);
	
	B = y[0][N - 1];
	s = s - (B - q_b)/Y;
	
	for (int i = 0; i < 2; i++){
		delete [] y[i];
	}
	delete [] y;
	
	while (abs(B - q_b) > eps and k <= 50){
		y = Runge_Kutta_4_order_full(f_nl, h, N, a, q_a, s);
		Y = RK_4_YU(2 * h, (N - 1)/2 + 1, a, 0, 1, y[0], y[1]);
		
		B = y[0][N - 1];
		s = s - (B - q_b)/Y;
		
		for (int i = 0; i < 2; i++){
			delete [] y[i];
		}
		delete [] y;
		
		k++;
	}
	y_f = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s);
	
	for (int i = 0; i < 2; i++){
		delete [] y[i];
	}
	delete [] y;
	
	return y_f;
}


double* Secant(double h, int N, double a, double q_a, double q_b){
	int k = 1;
	double s = 0.5, Y, eps = 1e-4, delta = 1e-4, B_1, B_2;
	double *y_1 = new double[N];
	double *y_2 = new double[N];
	
	
	y_1 = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s);
	y_2 = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s + delta);
	
	B_1 = y_1[N - 1];
	B_2 = y_2[N - 1];
	s = s - (B_1 - q_b) * delta / (B_2 - B_1);
	
	delete [] y_1;
	delete [] y_2;
	
	while (abs(B_1 - q_b) > eps and k <= 50){
		y_1 = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s);
		y_2 = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s + delta);
	
		B_1 = y_1[N - 1];
		B_2 = y_2[N - 1];
		s = s - (B_1 - q_b) * delta / (B_2 - B_1);
		
		delete [] y_1;
		delete [] y_2;
		
		k++;
	}
	y_1 = Runge_Kutta_4_order(f_nl, h, N, a, q_a, s);
	
	delete [] y_2;
	
	return y_1;
}


void error(double* Y_0, double* Y_1, int N){
	double L_1_y_abs = 0; 		// L_1 абсолютное отклонение h относительно h/2
	double L_2_y_abs = 0; 		// L_2 абсолютное отклонение h относительно h/2
	double L_i_y_abs = 0;		// L_i абсолютное отклонение h относительно h/2
	double L_1_y = 0;			// L_1 норма решения y на сетке с шагом h
	double L_2_y = 0;			// L_2 норма решения y на сетке с шагом h
	double L_i_y = 0;			// L_i норма решения y на сетке с шагом h
	
	//ofstream fout("errors.txt");
	//fout.setf(ios::scientific);
	
	
	for (int i = 0; i < N; i++){
		L_1_y 	  += abs(Y_0[i]);
		L_1_y_abs += abs(Y_0[i] - Y_1[2 * i]);
		
		L_2_y 	  += Y_0[i] * Y_0[i];
		L_2_y_abs += (Y_0[i] - Y_1[2 * i]) * (Y_0[i] - Y_1[2 * i]);
		
		L_i_y 	  = max(L_i_y, Y_0[i]);
		L_i_y_abs = max(L_i_y_abs, abs(Y_0[i] - Y_1[2 * i]));
	}
	
	L_2_y_abs = sqrt(L_2_y_abs);
	L_2_y 	  = sqrt(L_2_y);
	
	//cout << L_i_y << " " << L_i_z;
	
	cout.setf(ios::scientific);
	
	cout << "              |-----------------------------|" << endl;
	cout << "              |   Absolute   |   Relative   |" << endl;
	cout << "              |      h/2     |      h/2     |" << endl;
	cout << "L_1 error     |" << setw(14) << L_1_y_abs << '|' << setw(14) << L_1_y_abs / L_1_y << '|' << endl;
	cout << "L_2 error     |" << setw(14) << L_2_y_abs << '|' << setw(14) << L_2_y_abs / L_2_y << '|' << endl;
	cout << "L_i error     |" << setw(14) << L_i_y_abs << '|' << setw(14) << L_i_y_abs / L_i_y << '|' << endl;
	cout << "              |-----------------------------|" << endl;
	
}


int main(void){
	double a = 0., b = 1., q_a_l = -3.5, q_b_l = 3., q_a_nl = -0.5, q_b_nl = 0.5, h;
	int N = 55;
	double *Y_0;
	double *Y_1;
	
	h = (b - a) / (N - 1);
	
	ofstream fout1("1.txt");
	ofstream fout2("2.txt");
	ofstream fout3("3.txt");
	ofstream fout4("4.txt");
	fout1.setf(ios::scientific);
	fout2.setf(ios::scientific);
	fout3.setf(ios::scientific);
	fout4.setf(ios::scientific);
	cout.setf(ios::scientific);
	
	
	Y_0 = Shooting(h,   N,         a, q_a_l, q_b_l);
	Y_1 = Shooting(h/2, 2 * N - 1, a, q_a_l, q_b_l);
	for (int i = 0; i < N; i++){
		fout1 << a + h * i << " " << Y_0[i] << endl;
	}
	cout << "              |-----------------------------|" << endl;
	cout << "              |           Shooting          |" << endl;
	error(Y_0, Y_1, N);
	
	delete [] Y_0;
	delete [] Y_1;
	
	
	Y_0 = Tridiagonal(h,   N,         a, q_a_l, q_b_l);
	Y_1 = Tridiagonal(h/2, 2 * N - 1, a, q_a_l, q_b_l);
	for (int i = 0; i < N; i++){
		fout2 << a + h * i << " " << Y_0[i] << endl;
	}
	cout << "              |         Tridiagonal         |" << endl;
	error(Y_0, Y_1, N);
	
	delete [] Y_0;
	delete [] Y_1;
	
	
	Y_0 = Newton(h,   N,         a, q_a_nl, q_b_nl);
	Y_1 = Newton(h/2, 2 * N - 1, a, q_a_nl, q_b_nl);
	for (int i = 0; i < N; i++){
		fout3 << a + h * i << " " << Y_0[i] << endl;
	}
	cout << "              |            Newton           |" << endl;
	error(Y_0, Y_1, N);
	
	delete [] Y_0;
	delete [] Y_1;
	
	
	Y_0 = Secant(h,   N,         a, q_a_nl, q_b_nl);
	Y_1 = Secant(h/2, 2 * N - 1, a, q_a_nl, q_b_nl);
	for (int i = 0; i < N; i++){
		fout4 << a + h * i << " " << Y_0[i] << endl;
	}
	cout << "              |            Secant           |" << endl;
	error(Y_0, Y_1, N);
	
	delete [] Y_0;
	delete [] Y_1;
	
	
	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	
	system("py task4.py");
	
	return 0;
}
