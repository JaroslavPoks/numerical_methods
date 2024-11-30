#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;


double f(double x){
	//return sin(x);
	return x * x * x;
}


double I_f(double a, double b){
	return b * b * b * b / 4. - a * a * a * a / 4.;
}


double rect(double x_l, double x_r, double h){					// Формула прямоугольников
	return f((x_l + x_r) / 2.) * h;
}


double trap(double x_l, double x_r, double h){					// Формула трапеций
	return (f(x_l) + f(x_r)) * h / 2.;
}


double Simpson(double x_l, double x_r, double h){				// Формула Симпсона
	return (f(x_l) + 4. * f((x_l + x_r) / 2.) + f(x_r)) * h / 6.;
}	


double NC(double x_l, double x_r, double h){					// Формула Ньютона-Котеса
	return (19. * f(x_l) + 75. * f(x_l + (x_r - x_l) / 5.) + 50. * f(x_l + 2. * (x_r - x_l) / 5.) + 
			50. * f(x_l + 3. * (x_r - x_l) / 5.) + 75. * f(x_l + 4. * (x_r - x_l) / 5.) + 19. * f(x_r)) * h / 288.;
}


double Gauss(double x_l, double x_r, double h){					// Формула Гаусса
	return (5. * f((x_r + x_l) / 2. + (x_r - x_l) * (-sqrt(0.6)) / 2.) + 8. * f((x_r + x_l) / 2.) + 
			5. * f((x_r + x_l) / 2. + (x_r - x_l) * sqrt(0.6) / 2.)) * h / 18.;
}


int main(void){
	int K = 10, Z = 300;
	double a = 0.0, b = 1.0, I_rect = 0.0, I_trap = 0.0, I_NC = 0.0, I_Gauss = 0.0, I_Simpson = 0.0, h;
	double rel_err_rect_1 = 0.0, rel_err_trap_1 = 0.0, rel_err_NC_1 = 0.0, rel_err_Gauss_1 = 0.0, rel_err_Simpson_1 = 0.0;
	double rel_err_rect_2 = 0.0, rel_err_trap_2 = 0.0, rel_err_NC_2 = 0.0, rel_err_Gauss_2 = 0.0, rel_err_Simpson_2 = 0.0;
	double tmp = 0.0;
	
	h = (b - a) / K;
	ofstream fout("1.csv");
	cout.setf(ios::scientific);
	fout.setf(ios::scientific);
	
	for (int i = 0; i < K; i++){											// Вычисление интегралов и норм на K отрезках
		I_rect 	  += rect(a + i * h, a + (i + 1) * h, h);
		I_trap	  += trap(a + i * h, a + (i + 1) * h, h);
		I_Simpson += Simpson(a + i * h, a + (i + 1) * h, h);
		I_NC 	  += NC(a + i * h, a + (i + 1) * h, h);
		I_Gauss   += Gauss(a + i * h, a + (i + 1) * h, h);
	}
	
	tmp = I_f(a, b);
	rel_err_rect_1    = abs(I_rect - tmp) / abs(tmp);
	rel_err_trap_1    = abs(I_trap - tmp) / abs(tmp);
	rel_err_Simpson_1 = abs(I_Simpson - tmp) / abs(tmp);
	rel_err_NC_1      = abs(I_NC - tmp) / abs(tmp);
	rel_err_Gauss_1   = abs(I_Gauss - tmp) / abs(tmp);
	
	I_rect    = 0.;
	I_trap 	  = 0.;
	I_Simpson = 0.;
	I_NC	  = 0.;
	I_Gauss   = 0.;
	h = (b - a) / (2 * K);
	
	for (int i = 0; i < 2 * K; i++){										// Вычисление интегралов и норм на 2K отрезках
		I_rect 	  += rect(a + i * h, a + (i + 1) * h, h);
		I_trap	  += trap(a + i * h, a + (i + 1) * h, h);
		I_Simpson += Simpson(a + i * h, a + (i + 1) * h, h);
		I_NC 	  += NC(a + i * h, a + (i + 1) * h, h);
		I_Gauss   += Gauss(a + i * h, a + (i + 1) * h, h);
	}
	
	rel_err_rect_2    = abs(I_rect - tmp) / abs(tmp);
	rel_err_trap_2    = abs(I_trap - tmp) / abs(tmp);
	rel_err_Simpson_2 = abs(I_Simpson - tmp) / abs(tmp);
	rel_err_NC_2      = abs(I_NC - tmp) / abs(tmp);
	rel_err_Gauss_2   = abs(I_Gauss - tmp) / abs(tmp);
	
	fout << K << "\n";
	fout << rel_err_rect_1    << "," << rel_err_rect_2    << "\n";
	fout << rel_err_trap_1    << "," << rel_err_trap_2    << "\n";
	fout << rel_err_Simpson_1 << "," << rel_err_Simpson_2 << "\n";
	fout << rel_err_NC_1      << "," << rel_err_NC_2      << "\n";
	fout << rel_err_Gauss_1   << "," << rel_err_Gauss_2   << "\n";
	
	//cout << I_f(a, b) << " " << I_rect << " " << I_trap << " " << I_Simpson << " " << I_NC << " " << I_Gauss;
	fout.close();
	system("py task5.py");
	
	return 0;
}
	