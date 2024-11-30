#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;


double Phi_1(double y){
	return sqrt(y + 1);
}


double Phi_2(double x){
	return atan(x);
}


double f_1(double x, double y){
	return x * x - 1 - y;
}


double f_2(double x, double y){
	return y - atan(x);
}


double determinant(double x){
	return (2 * x + 2 * x * x * x - 1) / (1 + x * x);
}


int main(void){
	double x_0 = -1., y_0 = -1., r = 0., x, y, x_p, y_p, dx, dy, det, eps = 1e-5;
	int k = 0;
	
	
	x = Phi_1(y_0);																				// Метод простых итераций
	y = Phi_2(x_0);
	x_p = x;
	y_p = y;
	r = sqrt((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0));
	
	while (r > eps and k <= 500){
		x = Phi_1(y_p);
		y = Phi_2(x_p);
		r = sqrt((x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));
		x_p = x;
		y_p = y;
		k++;
	}
	
	cout << "Fixed-point iteration: x = " << x << ", y = " << y << ", Number of iterations: " << k << ", r = " << r << endl;
	
	x   = 0;
	y   = 0;
	k   = 0;
	r   = 0;
	x_p = 0;
	y_p = 0;
	
	det = determinant(x_0);
	dx  = - (f_1(x_0, y_0) + f_2(x_0, y_0)) / det;																		// Метод Ньютона
	dy  = - (f_1(x_0, y_0) / (1 + x_0 * x_0) + f_2(x_0, y_0) * 2 * x_0) / det;
	x   = x_0 + dx;
	y   = y_0 + dy;
	x_p = x;
	y_p = y;
	r   = sqrt((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0));
	k++;
	
	while  (r > eps and k <= 500){
		det = determinant(x_p);
		dx  = - (f_1(x_p, y_p) + f_2(x_p, y_p)) / det;
		dy  = - (f_1(x_p, y_p) / (1 + x_p * x_p) + f_2(x_p, y_p) * 2 * x_p) / det;
		r   = sqrt((x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));
		x   = x_p + dx;
		y   = y_p + dy;
		x_p = x;
		y_p = y;
		k++;
	}
	
	cout << "Newton's method: x = " << x << ", y = " << y << ", Number of iterations: " << k << ", r = " << r << endl;
	
	x   = 0;
	y   = 0;
	k   = 0;
	r   = 0;
	x_p = 0;
	y_p = 0;
	
	
	det = determinant(x_0);
	dx  = - (f_1(x_0, y_0) + f_2(x_0, y_0)) / det;																		// Модифицированный метод Ньютона
	dy  = - (f_1(x_0, y_0) / (1 + x_0 * x_0) + f_2(x_0, y_0) * 2 * x_0) / det;
	x   = x_0 + dx;
	y   = y_0 + dy;
	x_p = x;
	y_p = y;
	r   = sqrt((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0));
	k++;
	
	while  (r > eps and k <= 500){
		dx = - (f_1(x_p, y_p) + f_2(x_p, y_p)) / det;
		dy = - (f_1(x_p, y_p) / (1 + x_p * x_p) + f_2(x_p, y_p) * 2 * x_p) / det;
		r = sqrt((x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));
		x = x_p + dx;
		y = y_p + dy;
		x_p = x;
		y_p = y;
		k++;
	}
	
	cout << "Modified Newton's method: x = " << x << ", y = " << y << ", Number of iterations: " << k << ", r = " << r << endl;
	
	x   = 0;
	y   = 0;
	k   = 0;
	r   = 0;
	x_p = 0;
	y_p = 0;
	
	
	double df1_dx, df2_dx, df1_dy, df2_dy, h = 0.1;
	df1_dx = (f_1(x_0 + h, y_0) - f_1(x_0, y_0)) / h;
	df2_dx = (f_2(x_0 + h, y_0) - f_2(x_0, y_0)) / h;
	df1_dy = (f_1(x_0, y_0 + h) - f_1(x_0, y_0)) / h;
	df2_dy = (f_2(x_0, y_0 + h) - f_2(x_0, y_0)) / h;
	
	det = df1_dx * df2_dy - df2_dx * df1_dy;
	dx  = - ( f_1(x_0, y_0) * df2_dy - f_2(x_0, y_0) * df1_dy) / det;													// Дискретный метод Ньютона
	dy  = - (-f_1(x_0, y_0) * df2_dx + f_2(x_0, y_0) * df1_dx) / det;
	x   = x_0 + dx;
	y   = y_0 + dy;
	x_p = x;
	y_p = y;
	r   = sqrt((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0));
	k++;
	
	while  (r > eps and k <= 500){
		df1_dx = (f_1(x_p + h, y_p) - f_1(x_p - h, y_p)) / (2 * h);
		df2_dx = (f_2(x_p + h, y_p) - f_2(x_p - h, y_p)) / (2 * h);
		df1_dy = (f_1(x_p, y_p + h) - f_1(x_p, y_p - h)) / (2 * h);
		df2_dy = (f_2(x_p, y_p + h) - f_2(x_p, y_p - h)) / (2 * h);
	
		det = df1_dx * df2_dy - df2_dx * df1_dy;
		dx  = - ( f_1(x_p, y_p) * df2_dy - f_2(x_p, y_p) * df1_dy) / det;
		dy  = - (-f_1(x_p, y_p) * df2_dx + f_2(x_p, y_p) * df1_dx) / det;
		x   = x_p + dx;
		y   = y_p + dy;
		x_p = x;
		y_p = y;
		r   = sqrt((x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));
		k++;
	}
	
	cout << "Discrete Newton's method: x = " << x << ", y = " << y << ", Number of iterations: " << k << ", r = " << r << endl;
	
	system("py task2.py");
	
	return 0;
}
