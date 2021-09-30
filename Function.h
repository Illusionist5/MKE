#pragma once
#include "Struct.h"

//функция
double function(double x, double y, double t)
{
	return 0;
}

//коэфицент дифузии (лямбда) или коэфицент теплопроводности
const double dif_coef = 1;

//коэфицент плотности объемного телпа (гамма)
const double gamma_coef = 1;

//Краевое условие
double Kraev_1(int num_f, double x, double y, double t)
{
	switch (num_f)
	{
	case 0:
		return 1;
	case 1:
		return 0;
	default: 0;
	}
	return 0;
}

//Подсчет краевого условия первого рода
void uc_kr1()
{
	ifstream fin("Cond1.txt");
	int p, v1, v2;		//кол-во ребер, первая вершина, вторая вершина, номер краевой функции (0 по стандарту)
	double num_f;
	long B = 1e+10;
	double x1, y1, x2, y2;
	fin >> p;
	for (int i = 0; i < p; i++)
	{
		num_f = 0;
		fin >> v1;
		di[v1] = B;
		global_vector[v1] = B * Kraev_1(num_f, nod[v1].x, nod[v1].y, times[sloy]);
		/*
		x1 = nod[v1].x;
		y1 = nod[v1].y;
		x2 = nod[v2].x;
		y2 = nod[v2].y;
		double res1 = Kraev_1(num_f, x1, y1);
		double res2 = Kraev_1(num_f, x2, y2);
		//Записываем в матрицу
		di[v1] = B;
		//di[v2] = B;
		//Записываем в вектор
		global_vector[v1] = B * num_f;
		//global_vector[v2] = B * res2;
		*/
	}

	for (int i = 4; i < 8; i++)
	{
		global_vector[i] = 0;
		di[i] = B;
	}
		
}

/*
void read_cond(string cond_file_1, string cond_file_2) {
	ifstream fin1(cond_file_1 + ".txt");
	fin1 >> ncond1;
	cond1.resize(ncond1);
	for (int i = 0; i < ncond1; i++) {
		cond1[i].resize(4);
		for (int j = 0; j < 4; j++)
			fin1 >> cond1[i][j];
	}
	fin1.close();

	ifstream fin2(cond_file_2 + ".txt");
	fin2 >> ncond2;
	cond2.resize(ncond2);
	for (int i = 0; i < ncond2; i++) {
		cond2[i].resize(5);
		for (int j = 0; j < 5; j++)
			fin2 >> cond2[i][j];
	}
	fin2.close();
}
*/


// Ниже решение матрицы
//---------------------------------------------------//

vector<double> Ma;

void Mult(vector<double>& x, vector<double>& Ma)	//умножение на А
{
	for (int i = 0; i < n; ++i) 
	{
		int gi = ig[i], gi_1 = ig[i + 1];
		Ma[i] = di[i] * x[i];
		for (int j = gi; j < gi_1; ++j) 
		{
			int column = jg[j];
			Ma[i] += ggl[j] * x[column];
			Ma[column] += ggl[j] * x[i];
		}
	}
}

double Mult_scal(vector<double> vec1, vector<double> vec2)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

void LOS()
{
	int k = 0;
	double norm_eps;
	vector <double> r0(n, 0);
	vector <double> z0(n, 0);
	vector <double> p0(n, 0);
	vector <double> rk(n, 0);
	vector <double> zk(n, 0);
	vector <double> pk(n, 0);
	vector <double> xk(n, 0);
	Ma.resize(n);
	double maxiter = 100;
	double eps = 1e-12;

	for (int i = 0; i < n; ++i)
	{
		b[i] = 1;		//Ax0
	}
	int count = 0;
	Mult(b, Ma);
	for (int i = 0; i < n; ++i)
	{
		r0[i] = global_vector[i] - Ma[i];	//f - Ax0
		z0[i] = r0[i];
	}
	Mult(z0, p0);
	double sr = Mult_scal(r0, r0);
	while (sr > eps && count <= maxiter)
	{
		double pp = Mult_scal(p0, p0);
		double ak = Mult_scal(p0, r0) / pp;
		for (int i = 0; i < n; ++i)
		{
			b[i] = b[i] + ak * z0[i];
			r0[i] = r0[i] - ak * p0[i];
		}
		Mult(r0, Ma);
		double bk = -Mult_scal(p0, Ma) / pp;
		for (int i = 0; i < n; ++i)
		{
			z0[i] = r0[i] + bk * z0[i];
			p0[i] = Ma[i] + bk * p0[i];
		}
		sr = sqrt(Mult_scal(r0, r0));
		++count;

	}
	cout << endl;
	cout << "norm_eps = " << eps << endl;
	/*for (int i = 0; i < n; i++)
		cout << b[i] << '\n';
		*/
}
