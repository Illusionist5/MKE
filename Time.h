#pragma once
#include "Struct.h"

// ������ �����
void read_time()
{
	ifstream fin("Time.txt");	
	fin >> n_time;
	
	times.resize(n_time);
	fin >> times[0] >> times[n_time - 1];	// ��������� ������ � �������� �������, �.�. �������, ��� ��������� ����� ��������� ������ ����������
	fin.close();
	dt = (times[n_time - 1] - times[0]) / (n_time - 1);

	for (int i = 1; i < n_time; i++)
		times[i] = times[i - 1] + dt;
}

// ������ ������ ���� q1
void Q(int qn)
{
	for (int i = 0; i < n; i++)
		q[qn][i] = b[i];
}

// ������� M � ���������� ����������� �������
void M_Global(int el_num)	//������ ���������� ������� (� ������� �������, ������� ���������)
{
	int gl_ele[6];
	for (int i = 0; i < 6; i++)
		gl_ele[i] = el[el_num].n_nodes[i];

	for (int i = 0; i < 6; i++)
	{
		Rdi[gl_ele[i]] += lM[i][i];	//������������ �������
		for (int j = 0; j < i; j++)
		{
			auto a = gl_ele[i];
			auto b = gl_ele[j];
			if (a < b) swap(a, b);	//��� ����� � � �

			auto begin = jg.begin() + ig[a];
			if (ig[a + 1] > ig[a])
			{
				auto end = jg.begin() + ig[a + 1] - 1;
				auto iter = lower_bound(begin, end, b); //���������
				auto index = iter - jg.begin();
				Rggl[index] += lM[i][j];
			}
		}
	}
}

void Mult_M_global(vector<double>& x, vector<double>& Mm)	//��������� �� M
{
	for (int i = 0; i < n; ++i)
	{
		int gi = ig[i], gi_1 = ig[i + 1];
		Mm[i] = di[i] * x[i];
		for (int j = gi; j < gi_1; ++j)
		{
			int column = jg[j];
			Mm[i] += ggl[j] * x[column];
			Mm[column] += ggl[j] * x[i];
		}
	}
}

void Sub(vector<double>& x, vector<double>& y)	// �������� �� x y, ��������� � x
{
	for (int i = 0; i < x.size(); i++)
		x[i] -= y[i];
}

