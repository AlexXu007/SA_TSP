#include "HeuristicOperator.h"
#include<fstream>	//ifstream()��ofstream()
#include<vector>	//DM
#include<cmath>		//sqrt()��pow()
#include<cstdlib>	//rand()
#include<chrono>	//srand()
using namespace std;

//����1����ȡԭʼ�����Լ��������
vector<vector<double>> HeuristicOperator::getCoord(void) {
	//���ļ�
	ifstream originalData("data.txt");
	if (!originalData) {
		cerr << "File couldn't be opened" << endl;
		exit(EXIT_FAILURE);
	}

	double x;
	double y;
	vector<double> X;
	vector<double> Y;

	//��ȡ���ݲ���ֵ
	while (originalData >> x >> y) {
		X.push_back(x);
		Y.push_back(y);
	}
	vector<vector<double>> Coord(2, vector<double>(int(X.size())));
	Coord[0] = X;
	Coord[1] = Y;

	return Coord;
}

//����2������������
vector<vector<double>> HeuristicOperator::getDM(vector<vector<double>> Coord) {
	int n = int(Coord[0].size());
	vector<vector<double> > DM(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			DM[i][j] = sqrt(pow(Coord[0][i] - Coord[0][j], 2) + pow(Coord[1][i] - Coord[1][j], 2));
		}
	}
	return DM;
}

//����3�����ɳ�ʼ��
vector<int> HeuristicOperator::getInitS(int n) {
	vector<int> S0(n);
	//����˳������
	for (int i = 0; i < n; i++) {
		S0[i] = i;
	}
	//����˳������
	unsigned int seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	shuffle(S0.begin(), S0.end(), std::default_random_engine(seed));

	return S0;
}

//����4�����ۺ���
double HeuristicOperator::Eval(vector<int> S, vector<vector<double>> DM, int n) {
	double f = 0.0;
	for (int i = 0; i < n - 1; i++) {
		f = f + DM[S[i]][S[i + 1]];
	}
	f = f + DM[S[n - 1]][S[0]];

	return f;
}

//����5��������Χ����������ֵ������Ӧ��λ�ú���
vector<double> HeuristicOperator::bestS(vector<double> Eval, int Length) {
	vector<double> bS(2);		//�洢bestEval �� bestIndex
	double bestEval = 1e9;		//��������ֵ
	double bestIndex = -1;		//��������ֵ��Ӧ��λ��

	for (int i = 0; i < Length; i++) {
		if (Eval[i] < bestEval) {
			bestEval = Eval[i];
			bestIndex = (double)i;
		}
	}
	bS[0] = bestEval;
	bS[1] = bestIndex;

	return bS;
}

//����6������Sharking����λ�ú���
vector<int> HeuristicOperator::RandPosition(int n) {
	vector<int> RP(2);
	unsigned int seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	RP[0] = rand() % n;
	RP[1] = RP[0];

	while (!(RP[0] - RP[1])) {
		RP[1] = rand() % n;
	}

	return RP;
}

//����7����������
vector<int> HeuristicOperator::Swap(vector<int> S, vector<int> RP) {
	vector<int> newS = S;
	//��������
	int tmp = newS[RP[0]];
	newS[RP[0]] = newS[RP[1]];
	newS[RP[1]] = tmp;

	return newS;
}

//����8����ת����
vector<int> HeuristicOperator::Flip(vector<int> S, vector<int> RP) {
	vector<int> newS = S;
	//��ת����
	for (int i = 0; i<int(S.size()); i++) {
		if (i >= RP[0] && i <= RP[1]) {
			newS[i] = S[RP[0] + RP[1] - i];
		}
	}

	return newS;
}

//����9����������
vector<int> HeuristicOperator::Insert(vector<int> S, vector<int> RP) {
	vector<int> newS = S;
	//�������
	int tmp = newS[RP[0]];
	newS.erase(newS.begin() + RP[0]);
	newS.insert(newS.begin() + RP[1], tmp);

	return newS;
}