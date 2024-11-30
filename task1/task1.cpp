#include <iostream>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <fstream>

using namespace std;


double f(double x){
	return sin(x);
	//return 0.15 *  x * x * x * x + 2 * x * x + 1;
}


double Lagrange(double* points, double* denom, int size, double x){				// Многочлен Лагранжа
	double value_x = 0.0, num = 1.0;
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			if (j != i){
				num *= x - points[j];
			}
		}
		value_x += f(points[i]) * num / (denom[i]);
		num = 1.0;
	}
	return value_x;
}


int main(void){
	double a = -1.0, b = 5.0, delta, epsilon, den = 1.0;
	int N = 16, M = 200, k = 1;
	double* pointsN;
	double* denoms;
	
	denoms = new double[N];
	pointsN = new double[N];
	delta = (b - a) / (M - 1);
	epsilon = 0.5 * (b - a) / (N - 1);
	
	ofstream fout("1.txt");
	fout << a << " " << b << " " << N << " " << M << endl;
	
	uniform_real_distribution<double> distribution{a, b};						// Генератор случайных чисел
	mt19937 generator{random_device().operator ()()};
	
	pointsN[0] = a;
	while (k < N){																// Заполнение массива
		pointsN[k] = distribution(generator);
		for (int j = 0; j < k - 1; j++){
			if (abs(pointsN[k] - pointsN[j]) < epsilon)
				k--;
		}
		k++;
	}
	sort(pointsN, pointsN + N);
	pointsN[0] = a;
	pointsN[N - 1] = b;
	
	
	for (int i = 0; i < 2 * N; i++){
		if (i < N)
			fout << f(pointsN[i]) << " ";
		else if (i == N){
			fout << endl << pointsN[0] << " ";
		}
		else
			fout << pointsN[i - N] << " ";
	}
	fout << endl;
	
	for (int i = 0; i < N; i++){												// Вычисление знаменателей
		for (int j = 0; j < N; j++){
			if (j != i)
				den *= (pointsN[i] - pointsN[j]);
		}
		denoms[i] = den;
		den = 1.0;
	}
	

	for (int i = 0; i < M; i++){												// Построение многочлена Лагранжа
		fout << Lagrange(pointsN, denoms, N, a + i * delta) << " ";
	}
	
	fout.close();
	delete [] pointsN;
	delete [] denoms;
	
	system("py task1.py");
	
	return 0;
}
	