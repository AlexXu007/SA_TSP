#pragma once
#include<iostream>
#include<vector>
#include <algorithm>		// std::shuffle
#include <random>			// std::default_random_engine
#include<chrono>
using namespace std;

class HeuristicOperator {
public:
	vector<vector<double>> getCoord(void);		//����1����ȡ���꺯��
	vector<vector<double>> getDM(vector<vector<double>> Coord);		//����2����ȡ���������
	vector<int> getInitS(int n);		//����3����ȡ��ʼ�⺯��
	double Eval(vector<int> S, vector<vector<double>> DM, int n);	//����4�����ۺ���

	vector<double> bestS(vector<double> Eval, int Length);	//����5��������Χ����������ֵ������Ӧ��λ�ú���

	vector<int> RandPosition(int n);		//����6������Sharking����λ�ú���
	vector<int> Swap(vector<int> S, vector<int> RP);	//����7����������
	vector<int> Flip(vector<int> S, vector<int> RP);	//����8����ת����
	vector<int> Insert(vector<int> S, vector<int> RP);	//����9����������
};