//#include "Struct.h"
#include "Function.h"
#include "Time.h"

void Input()
{
	ifstream fin("input.txt");	//��������� ����, � ������� ���������� ������
	fin >> n >> m;	//����, ����� � ���, ��� n ���������

	//������������ �������
	nod.resize(n);
	el.resize(m);

	//������ ����
	for (int i = 0; i < n; i++)
		fin >> nod[i].x >> nod[i].y;

	//������ �������
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 6; j++)		//������ ������ ������� �����
			fin >> el[i].n_nodes[j];	//������ ����� ���������� � 0, ��� � � ���������� ����

		/*// �������������� ��������� �������
		nodes input[3];
		input[0].x = (nod[el[i].n_nodes[0]].x + nod[el[i].n_nodes[0]].x) / 2;
		input[0].y = (nod[el[i].n_nodes[0]].y + nod[el[i].n_nodes[0]].y) / 2;

		input[1].x = (nod[el[i].n_nodes[1]].x + nod[el[i].n_nodes[2]].x) / 2;
		input[1].y = (nod[el[i].n_nodes[1]].y + nod[el[i].n_nodes[2]].y) / 2;

		input[2].x = (nod[el[i].n_nodes[2]].x + nod[el[i].n_nodes[0]].x) / 2;
		input[2].y = (nod[el[i].n_nodes[2]].y + nod[el[i].n_nodes[0]].y) / 2;

		for (int j = 0; j < 3; j++)		//����������� �������� ��������
		{
			nod.push_back(input[i]);		//��������� � ����� ������ �����
			el[i].n_nodes[j+3] = nod.size() - 1;	//����������� �����
		}
		*/
		fin >> el[i].number;	//������� ���������� � ����
	}
	//n = nod.size();		//������ ���-�� ����� �� �������

	//��������� ����� � ������ times[n_time]
	read_time();
		
	fin.close();
}

void Gen_portrait()	//��������� ig � jg
{
	vector <set<int>> list(n);	//������ �������

	//��������� ���� ������
	for (int ielem = 0; ielem < m; ielem++)
		for (int i = 0; i < 6; i++) //6 - ���-�� �����, �������������� ��������� ��������
			for (int j = i + 1; j < 6; j++)
			{
				int ind1 = el[ielem].n_nodes[i];
				int ind2 = el[ielem].n_nodes[j];
				if (ind1 < ind2) swap(ind1, ind2);
				list[ind1].insert(ind2);
			}	

	//�������� �������� �� ������
	ig.resize(n + 1);
	ig[0] = ig[1] = 0;

	for (int i = 2; i < n + 1; i++) 
		ig[i] = ig[i - 1] + list[i - 1].size();

	jg.resize(ig[n]);
	for (int i = 1, k = 0; i < n; i++)
		for (int j : list[i]) 
		{
			jg[k] = j;
			k++;
		}

}
//���������� ��� ��������� �������� ���������� ��� ������� ��������� L1 � L2. ���� -1 - L ���
double IntegL(int l1, int l2)	
{
	if ((l1 == -1) || (l2 == -1))
		if (l1 == l2)
			return 2;	//��� �������
		else 
			return 6;	//��� ������ �������� l

	if (l1 == l2)
		return 12;	//��� l^2

	if (l1 != l2)
		return 24;	// ��� ������ l1 � l2

	return 0;
}


void Prot_Psi(int psi)
{
	psi++;
	switch (psi)
	{
	case 4:
	{
		psi1 = 0;
		psi2 = 1;
		break;
	}
	case 5:
	{
		psi1 = 0;
		psi2 = 2;
		break;
	}
	case 6:
	{
		psi1 = 1;
		psi2 = 2;
		break;
	}
	default:
	{
		psi1 = -1; 
		psi2 = -1;
		break;
	}
	}
}

//���������� �����, � ������� ��� ������������� L
int Prot_Line(int vert)
{
	switch (vert)
	{
	case 1:
		return 6;
	case 2:
		return 5;
	case 3:
		return 4;
	default:
		return 0;
	}
	return 0;
}

void Solution_local_matrix(int el_num)	//������ ��������� ������� (� ������� ����� ������� � �������)
{
	//��� �������� ������� ��������� x � y ����� � ������ local_x � local_y
	vector <double> local_x(6);
	vector <double> local_y(6);

	for (int i = 0; i < 6; i++)
	{
		local_x[i] = nod[el[el_num].n_nodes[i]].x;
		local_y[i] = nod[el[el_num].n_nodes[i]].y;
	}

	//double mult = dif * det / 2;

	//������� ��������� alfa
	vector<vector<double>> a(3);

	a[0].push_back(local_x[1] * local_y[2] - local_x[2] * local_y[1]);
	a[0].push_back(local_y[1] - local_y[2]);
	a[0].push_back(local_x[2] - local_x[1]);

	a[1].push_back(local_x[2] * local_y[0] - local_x[0] * local_y[2]);
	a[1].push_back(local_y[2] - local_y[0]);
	a[1].push_back(local_x[0] - local_x[2]);

	a[2].push_back(local_x[0] * local_y[1] - local_x[1] * local_y[0]);
	a[2].push_back(local_y[0] - local_y[1]);
	a[2].push_back(local_x[1] - local_x[0]);

	//DetD
	double det = abs((local_x[1] - local_x[0]) * (local_y[2] - local_y[0]) - (local_x[2] - local_x[0]) * (local_y[1] - local_y[0]));

	//������� ��������� G
	vector <vector <double>> G;	//��������� ������� ���������
	//������� ��
	G.resize(6);
	for (int i = 0; i < 6; i++)
		G[i].resize(6);

	int i1, i2, j1, j2 = 0;	//��� ���������� ��� 4, 5 � 6 �� ��� ������ L��
	double sumA = 0;	//��� �������� ������, ����� ����
	//��� 1, 2 � 3
	for (int i = 0; i < 3; i++)
	{
		//� 1, 2 � 3
		for (int j = 0; j <= i; j++)
		{
			sumA = (a[i][1] * a[j][1]) + (a[i][2] * a[j][2]);
			G[i][j] = dif_coef * det * sumA * 0.5;
			G[j][i] = G[i][j];		//��������������
		}

		sumA = 0;
		// c 4, 5 � 6
		for (int j = 3; j < 6; j++)
		{
			Prot_Psi(j);	//������������ ��� �� ��� ���
			j1 = psi1; j2 = psi2;
			sumA = (a[i][1] * (a[j1][1] + a[j2][1])) + (a[i][2] * (a[j1][2] + a[j2][2]));
			G[i][j] = dif_coef * det * sumA;
			G[j][i] = G[i][j];
		}
	}

	sumA = 0;
	for (int i = 3; i < 6; i++)
		for (int j = 3; j <= i; j++)
		{
			//���������� ��� ���
			Prot_Psi(i);
			i1 = psi1; i2 = psi2;
			Prot_Psi(j);
			j1 = psi1; i2 = psi2;

			//sumA = (alpha[i1][1] * (alpha[j2][1] / IntegL(i2, j1) + alpha[j1][1] / IntegL(i2, j2))
			sumA = (a[i1][1] * a[j1][1] + a[i1][2] * a[j1][2]) / IntegL(i2, j2);
			sumA += (a[i1][1] * a[j2][1] + a[i1][2] * a[j2][2]) / IntegL(i2, j1);
			sumA += (a[i2][1] * a[j1][1] + a[i2][2] * a[j1][2]) / IntegL(i1, j2);
			sumA += (a[i2][1] * a[j2][1] + a[i2][2] * a[j2][2]) / IntegL(i1, j1);
			G[i][j] = dif_coef * det * sumA;
			G[j][i] = G[i][j];
		}
		/*
		double k1;
		double k2;
		double sumA1, sumA2;
		//c 4
		k1 = 0;
		k2 = 0;
		sumA1 = (alpha[i][1] * alpha[0][1]) + (alpha[i][2] * alpha[0][2]);
		sumA2 = (alpha[i][1] * alpha[1][1]) + (alpha[i][2] * alpha[1][2]);
		if (i == 1)
			k1 = 2 / 3;
		if (i == 0)
			k2 = 2 / 3;
		G[i][3] = dif_coef * det * (sumA1 * k1 + sumA2 *k2);
		G[3][i] = G[i][3];

		//� 5
		k1 = 0;
		k2 = 0;
		sumA1 = (alpha[i][1] * alpha[1][1]) + (alpha[i][2] * alpha[1][2]);
		sumA2 = (alpha[i][1] * alpha[2][1]) + (alpha[i][2] * alpha[2][2]);
		if (i == 1)
			k1 = 2 / 3;
		if (i == 2)
			k2 = 2 / 3; 
		G[i][4] = dif_coef * det * (sumA1 * k1 + sumA2 * k2);
		G[4][i] = G[4][i];

		//� 6
		k1 = 0;
		k2 = 0;
		sumA1 = (alpha[i][1] * alpha[0][1]) + (alpha[i][2] * alpha[0][2]);
		sumA2 = (alpha[i][1] * alpha[2][1]) + (alpha[i][2] * alpha[2][2]);
		if (i == 0)
			k1 = 2 / 3;
		if (i == 2)
			k2 = 2 / 3;
		G[i][5] = dif_coef * det * (sumA1 * k1 + sumA2 * k2);
		G[5][i] = G[5][i];
	}

	//������� 4 4
	sumA = pow(alpha[0][1], 2) + pow(alpha[0][2], 2) + pow(alpha[1][1], 2) + pow(alpha[1][2], 2);
	sumA = sumA + alpha[0][1] * alpha[1][1] + alpha[0][2] * alpha[1][2];
	G[3][3] = dif_coef * det * sumA * 4 / 3;

	//������� 5 5
	sumA = pow(alpha[1][1], 2) + pow(alpha[1][2], 2) + pow(alpha[2][1], 2) + pow(alpha[2][2], 2);
	sumA = sumA + alpha[1][1] * alpha[2][1] + alpha[1][2] * alpha[2][2];
	G[4][4] = dif_coef * det * sumA * 4 / 3;

	//������� 6 6
	sumA = pow(alpha[0][1], 2) + pow(alpha[0][2], 2) + pow(alpha[2][1], 2) + pow(alpha[2][2], 2);
	sumA = sumA + alpha[0][1] * alpha[2][1] + alpha[0][2] * alpha[2][2];
	G[5][5] = dif_coef * det * sumA * 4 / 3;

	//��� 4, 5 � 6
	for (int i = 3; i < 6; i++)
		for (int j = 3; j < i; j++)
		{
			G[i][j] = 1;
			G[j][i] = G[i][j];		//��������������m
		}
	*/
	cout << endl << "Matrix " << endl;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			cout << G[i][j] << " ";
		cout << endl;
	}

	//cout << det << endl;

	//������� ���� M
	vector <vector <double>> M;	//��������� ������� ����
	//������� ��
	M.resize(6);
	for (int i = 0; i < 6; i++)
		M[i].resize(6);

	//��� 1, 2 � 3
	for (int i = 0; i < 3; i++)
	{
		M[i][i] = gamma_coef * det / 12;	//������� ���������

		//c 1, 2 � 3
		for (int j = 0; j < i; j++)
		{
			M[i][j] = gamma_coef * det / 24;
			M[j][i] = M[i][j];		//��������������
		}

		//c 4, 5 � 6
		for (int j = 3; j < 6; j++)
		{
			M[i][j] = gamma_coef * det / 60;
			if (i == Prot_Line(i + 1))
				M[i][j] = M[i][j] / 2;
			M[j][i] = M[i][j];		//��������������
		}
	}
	//��� 4, 5 � 6
	for (int i = 3; i < 6; i++)
	{
		M[i][i] = gamma_coef * det / 180;	//���������
		for (int j = 3; j < i; j++)
		{
			M[i][j] = gamma_coef * det / 360;
			M[j][i] = M[i][j];	//��������������
		}
	}

	//�����!
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
		{
			M[i][j] = M[i][j] / dt;	// ������� �������� M � ������ �������
			lM[i][j] = M[i][j];	// ������ M ����������, ��� �������� ������ ����� � �����
		}
			
	
	/*
	cout << endl << "Matrix mass" << endl;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			cout << M[i][j] << " ";
		cout << endl;
	}
	*/

	//������� ��������� �������
	//������� ��
	local_matrix.resize(6);
	for (int i = 0; i < 6; i++)
		local_matrix[i].resize(6);

	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
		{
			local_matrix[i][j] = M[i][j] + G[i][j];	//�� �������
			if (sloy == 0)
				local_matrix[i][j] = local_matrix[i][j] * dt;
		}
			

	/*cout << endl <<"Local matrix" << endl;
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			cout << local_matrix[i][j] << " ";
		cout << endl;
	}
	*/

	//������� ������� ������ �����
	local_vector.resize(6);
	for (int i = 0; i < 6; i++)
	{
		local_vector[i] = 0;
		for (int j = 0; j < 6; j++)
		{
			//�� �������� ������ �� ��������� �����!
			local_vector[i] += function(local_x[j], local_y[j], 1) * M[i][j] / gamma_coef;	
		}
	}
		
}

void Global_matrix(int el_num)	//������ ���������� ������� (� ������� �������, ������� ���������)
{
	int gl_ele[6];
	for (int i = 0; i < 6; i++)
		gl_ele[i] = el[el_num].n_nodes[i];

	for (int i = 0; i < 6; i++)
	{
		di[gl_ele[i]] += local_matrix[i][i];	//������������ �������
		global_vector[gl_ele[i]] += local_vector[i];	//������ ������
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
				ggl[index] += local_matrix[i][j];
			}
		}
	}
}


int main()
{
	Input();
	Gen_portrait();
	
	di.resize(n);
	global_vector.resize(n);
	ggl.resize(ig[n] - ig[0]);

	//��� ����������� ������� M
	Rdi.resize(n);
	Rggl.resize(ig[n] - ig[0]);
	b.resize(n);

	//��� q
	q[0].resize(n);
	q[1].resize(n);
	q[2].resize(n);

	//��������� ������� M
	lM.resize(6);
	for (int i = 0; i < 6; i++)
		lM[i].resize(6);

	//Solution_local_matrix(0);

	sloy = 0;

	//������� ���������� �������
	for (int i = 0; i < m; i++)
	{
		Solution_local_matrix(i);
		Global_matrix(i);
		M_Global(i);
	}
	uc_kr1();

	LOS();
	Q(sloy);

	sloy++;
	
	vector<double> Mm;	//��������� ��������� ������� �� ���������� M
	Mm.resize(n);

	for (int i = 0; i < m; i++)
	{
		Solution_local_matrix(i);
		Global_matrix(i);
	}
	//uc_kr1();

	Mult_M_global(q[sloy], Mm);	//�����������, ��� �������� �� bi, ������ � MM
	Sub(global_vector, Mm);	//�������� �� bi

	uc_kr1();
	LOS();

	Q(sloy);
	Mm.clear();
	
	for (int j = 0; j < 9; j++)
		cout << q[sloy][j] << endl;

	cout << endl;

	/* 
	for (int i = 0; i < n + 1; i++)
	{
		cout << ig[i] << " " << jg[i];
		cout << endl;
	}
	*/

	cout << "All done\n" << endl;	//����� � ���������� �������
	system("pause");	//�������� ������� �������
	return 0;
}
