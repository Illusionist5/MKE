#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

struct nodes	//структура узлов
{
	double x;
	double y;
};

struct elems	//структура конечный элементов или областей
{
	int n_nodes[6];	//массив узлов приндалежащий конечному элементу
				//0,1,2 - угловые, 3,4,5 - середина
	int number;		//номер подобласти расчетной области
};

//глобальные переменные
int n;	//кол-во узлов
int m;	//кол-во областей или конечномерных элементов
vector <nodes> nod;	//массив узлов
vector <elems> el;	//массив конечномерных элементов

int psi1, psi2;	// для разложение пси 4, 5 и 6 на L

//для подсчета разряженной матрицы
vector <double> ggl, di;	//верхний треугольник и главная диогональ
vector <int> ig, jg;	//координаты, с которых начинаются ненулевые элементы

//локальные мемы
vector <vector <double>> local_matrix;	//локальная матрица
vector <double> local_vector;	//локальная матрица

//глобальные мемы
vector <vector <double>> global_A;
vector <double> global_vector;

//Время
int n_time;		//кол-во времен
vector <double> times;	//временная сетка
double dt;		

//Слои
vector <double> q[3];	//слои для вектора
int sloy;	//номер слоя, на котором находимся

//Вектор B
vector <double> b(n, 0);

// для расчетов
vector <vector <double>> lM;	// локальная матрица M
//для подсчета правой разряженной матрицы
vector <double> Rggl, Rdi;	//верхний треугольник и главная диогональ