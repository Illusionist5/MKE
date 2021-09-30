#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

struct nodes	//��������� �����
{
	double x;
	double y;
};

struct elems	//��������� �������� ��������� ��� ��������
{
	int n_nodes[6];	//������ ����� ������������� ��������� ��������
				//0,1,2 - �������, 3,4,5 - ��������
	int number;		//����� ���������� ��������� �������
};

//���������� ����������
int n;	//���-�� �����
int m;	//���-�� �������� ��� ������������� ���������
vector <nodes> nod;	//������ �����
vector <elems> el;	//������ ������������� ���������

int psi1, psi2;	// ��� ���������� ��� 4, 5 � 6 �� L

//��� �������� ����������� �������
vector <double> ggl, di;	//������� ����������� � ������� ���������
vector <int> ig, jg;	//����������, � ������� ���������� ��������� ��������

//��������� ����
vector <vector <double>> local_matrix;	//��������� �������
vector <double> local_vector;	//��������� �������

//���������� ����
vector <vector <double>> global_A;
vector <double> global_vector;

//�����
int n_time;		//���-�� ������
vector <double> times;	//��������� �����
double dt;		

//����
vector <double> q[3];	//���� ��� �������
int sloy;	//����� ����, �� ������� ���������

//������ B
vector <double> b(n, 0);

// ��� ��������
vector <vector <double>> lM;	// ��������� ������� M
//��� �������� ������ ����������� �������
vector <double> Rggl, Rdi;	//������� ����������� � ������� ���������