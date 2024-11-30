#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;


double f(double x){
	return sin(x);
	//return x * x * x;
}


double D_f(double x){
	return cos(x);
	//return 3 * x * x;
}


double D_central(double x, double h){														// Центральная производная
	return (f(x + h) - f(x - h)) / (2 * h);
}


double D_left(double x, double h){															// Левая производная
	return ( - 3 * f(x) + 4 * f(x + h) - f(x + 2 * h)) / (2 * h);
}


double D_right(double x, double h){															// Правая производная
	return (f(x - 2 * h) - 4 * f(x - h) + 3 * f(x)) / (2 * h);
}


int main(void){
	int M = 10, Z = 400;
	double h1, h2, a = -3.14, b = 3.14, R = 0.0, D_1 = 0.0, L_1_f = 0.0, L_2_f = 0.0, L_inf_f = 0.0;
	double L_1_rel_h  = 0.0, L_2_rel_h  = 0.0, L_inf_rel_h  = 0.0, L_1_abs_h  = 0.0, L_2_abs_h  = 0.0, L_inf_abs_h  = 0.0;
	double L_1_rel_h2 = 0.0, L_2_rel_h2 = 0.0, L_inf_rel_h2 = 0.0, L_1_abs_h2 = 0.0, L_2_abs_h2 = 0.0, L_inf_abs_h2 = 0.0;
	double L_1_rel_r  = 0.0, L_2_rel_r  = 0.0, L_inf_rel_r  = 0.0, L_1_abs_r  = 0.0, L_2_abs_r  = 0.0, L_inf_abs_r  = 0.0;
	
	h1 = (b - a) / (M - 1);
	h2 = (b - a) / (2 * (M - 1));
	
	ofstream fout1("1.txt");
	ofstream fout2("2.csv");	
	
	fout1 << a << " " << b << " " << Z << " " << M << endl;
	
	fout1.setf(ios::scientific);
	fout2.setf(ios::scientific);
	
	double tmp = 0.0;									// Вычисление производных и норм на сетке с шагом h/2
	D_1 = D_left(a, h2);
	tmp = D_f(a);
	fout1 << setw(13) << D_1 << endl;
	
	L_1_abs_h2 += abs(tmp - D_1);
	L_1_f += abs(tmp);
	
	L_2_abs_h2 += (tmp - D_1) * (tmp - D_1);
	L_2_f += tmp * tmp;
	
	L_inf_abs_h2 = abs(tmp - D_1);
	L_inf_f = abs(tmp);
	
	for (int i = 1; i < 2 * M - 2; i++){
		tmp = D_f(a + i * h2);
		D_1 = D_central(a + i * h2, h2);
		fout1 << setw(13) << D_1 << endl;
		
		L_1_abs_h2 += abs(tmp - D_1);
		L_1_f += abs(tmp);
		
		L_2_abs_h2 += (tmp - D_1) * (tmp - D_1);
		L_2_f += tmp * tmp;
		
		L_inf_abs_h2 = max(abs(tmp - D_1), L_inf_abs_h2);
		L_inf_f = abs(tmp);
	}
	tmp = D_f(b);
	D_1 = D_right(b, h2);
	fout1 << setw(13) << D_1 << endl;
	
	L_1_abs_h2 += abs(tmp - D_1);
	L_1_f += abs(tmp);
	
	L_2_abs_h2 += (tmp - D_1) * (tmp - D_1);
	L_2_f += tmp * tmp;
	
	L_inf_abs_h2 = max(abs(tmp - D_1), L_inf_abs_h2);
	L_inf_f = max(abs(tmp), L_inf_f);
	
	L_2_f = sqrt(L_2_f);
	L_2_abs_h2 = sqrt(L_2_abs_h2);
	L_2_rel_h2 = L_2_abs_h2 / L_2_f;
	
	L_1_rel_h2 = L_1_abs_h2 / L_1_f;
	
	L_inf_rel_h2 = L_inf_abs_h2 / L_inf_f;
	
	L_1_f = 0.0;
	L_2_f = 0.0;
	L_inf_f = 0.0;
	
	//////////////////////////////////////////////////////////
	
	D_1 = D_left(a, h1);													// Вычисление производных и норм и уточнение по Рунге на сетке с шагом h
	R = (D_1 - D_left(a, h2))/(pow(0.5, 2) - 1);
	tmp = D_f(a);
	fout1 << setw(13) << D_1 << endl;
	fout1 << setw(13) << R   << endl;
	
	L_1_abs_h += abs(tmp - D_1);
	L_1_abs_r += abs(tmp - D_1 - R);
	L_1_f     += abs(tmp);
	
	L_2_abs_h += (tmp - D_1) * (tmp - D_1);
	L_2_abs_r += (tmp - D_1 - R) * (tmp - D_1 - R);
	L_2_f     += tmp * tmp;
	
	L_inf_abs_h = abs(tmp - D_1);
	L_inf_abs_r = abs(tmp - D_1 - R);
	L_inf_f     = max(L_inf_f, abs(tmp));
	
	for (int i = 1; i < M - 1; i++){
		tmp = D_f(a + i * h1);
		D_1 = D_central(a + i * h1, h1);
		R = (D_1 - D_central(a + i * h1, h2))/(pow(0.5, 2) - 1);
		fout1 << setw(13) << D_1 << endl;
		fout1 << setw(13) << R   << endl;
		
		L_1_abs_h += abs(tmp - D_1);
		L_1_abs_r += abs(tmp - D_1 - R);
		L_1_f     += abs(tmp);
	
		L_2_abs_h += (tmp - D_1) * (tmp - D_1);
		L_2_abs_r += (tmp - D_1 - R) * (tmp - D_1 - R);
		L_2_f     += tmp * tmp;
	
		L_inf_abs_h = max(L_inf_abs_h, abs(tmp - D_1));
		L_inf_abs_r = max(L_inf_abs_r, abs(tmp - D_1 - R));
		L_inf_f     = max(L_inf_f, abs(tmp));
	}
	tmp = D_f(b);
	D_1 = D_right(b, h1);
	R = (D_1 - D_right(b, h2))/(pow(0.5, 2) - 1);
	fout1 << setw(13) << D_1 << endl;
	fout1 << setw(13) << R   << endl;
	
	L_1_abs_h += abs(tmp - D_1);
	L_1_abs_r += abs(tmp - D_1 - R);
	L_1_f     += abs(tmp);
	
	L_2_abs_h += (tmp - D_1) * (tmp - D_1);
	L_2_abs_r += (tmp - D_1 - R) * (tmp - D_1 - R);
	L_2_f     += tmp * tmp;
	
	L_inf_abs_h = max(L_inf_abs_h, abs(tmp - D_1));
	L_inf_abs_r = max(L_inf_abs_r, abs(tmp - D_1 - R));
	L_inf_f     = max(L_inf_f, abs(tmp));
	
	L_2_f = sqrt(L_2_f);
	L_2_abs_h = sqrt(L_2_abs_h);
	L_2_abs_r = sqrt(L_2_abs_r);
	L_2_rel_h = L_2_abs_h / L_2_f;
	L_2_rel_r = L_2_abs_r / L_2_f;
	
	L_1_rel_h = L_1_abs_h / L_1_f;
	L_1_rel_r = L_1_abs_r / L_1_f;
	
	L_inf_rel_h = L_inf_abs_h / L_inf_f;
	L_inf_rel_r = L_inf_abs_r / L_inf_f;
	
	
	fout2 << L_1_rel_h   << "," << L_1_rel_h2   << "," << L_1_rel_r   << "\n";
	fout2 << L_2_rel_h   << "," << L_2_rel_h2   << "," << L_2_rel_r   << "\n";
	fout2 << L_inf_rel_h << "," << L_inf_rel_h2 << "," << L_inf_rel_r << "\n";
	
	fout2 << L_1_abs_h   << "," << L_1_abs_h2   << "," << L_1_abs_r   << "\n";
	fout2 << L_2_abs_h   << "," << L_2_abs_h2   << "," << L_2_abs_r   << "\n";
	fout2 << L_inf_abs_h << "," << L_inf_abs_h2 << "," << L_inf_abs_r << "\n";
	
	fout1.close();
	fout2.close();
	
	system("py task4.py");
	
	return 0;
}