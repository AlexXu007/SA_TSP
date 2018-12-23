/*
	文件名：SA_TSP
	作者：Alex Xu
	地址：Dalian Maritime University
	描述：利用模拟退火算法求解TSP问题
	创建时间：2018年12月22日16:40
*/

#include<iostream>
#include<chrono>
#include<vector>
#include<algorithm>		//min()  max()
#include "HeuristicOperator.h"

using namespace std;
using namespace chrono;

int main() {

	//求解前准备工作
	auto start = system_clock::now();	//开始计时
	HeuristicOperator HO;	//引入各启发式算子
	vector<vector<double>> DM;
	DM = HO.getDM(HO.getCoord());	//生成距离矩阵
	int N = int(DM[0].size());	//城市规模

		//设置算法参数
	//double T0 ;
	double Tend = 0.01;
	double Alpha = 0.9;
	int MarkovL = 10*N;

	//初始化算法
	vector<int> S0(N);		//初始解
	vector<int> NS(N);		//当前解
	vector<int> BS(N);		//最优解
	double ns;		//记录当前代最优解的评价值
	vector<vector<int>> MK(MarkovL, vector<int>(N));	//记录Markov搜索过程中的解
	vector<double> mks(MarkovL);	//记录Markov搜索过程中解的评价值

	double bs = 1e9;				//记录全局最优解的评价值
	unsigned int seed = (unsigned)system_clock::now().time_since_epoch().count();	//产生随机种子

	//设置初始温度
	for (int i = 0; i < MarkovL; i++) {
		mks[i] = HO.Eval(HO.getInitS(N), DM, N);
	}

	auto MinMax = minmax_element(mks.begin(), mks.end());
	double T0 = (*MinMax.first - *MinMax.second) / log(0.9);
	cout << T0 << endl;
	double T = T0;			//当前温度
	//生成初始解并初始化当前解
	S0 = HO.getInitS(N);
	NS = S0;
	//计算当前解的评价值
	ns = HO.Eval(S0, DM, N);

	//开始模拟退火过程
	while (T > Tend) {

		//最优解的更新
		if (ns < bs) {
			bs = ns;
			BS = NS;
		} else{
			double P = rand() / double(RAND_MAX);
			double Pr = exp(-T/T0);

			cout << "当前温度为" << T
				<< "\t当前解为" << bs
				<< "\t此时Pr为" << Pr
				<< endl;

			if (P > Pr) {
				bs = ns;
				BS = NS;
			}
		}

		//解的扰动，Markov搜索
		for (int L = 0; L < MarkovL; L++) {
			vector<int> randPosition = HO.RandPosition(N);
			double randShk = rand() / double(RAND_MAX);
			if (randShk < 0.33) {
				MK[L] = HO.Swap(NS, randPosition);		//交换操作
				mks[L] = HO.Eval(MK[L], DM, N);			//记录该解的评价值
			}
			else if (randShk > 0.67) {
				MK[L] = HO.Flip(NS, randPosition);		//翻转操作
				mks[L] = HO.Eval(MK[L], DM, N);			//记录该解的评价值
			}
			else {
				MK[L] = HO.Insert(NS, randPosition);	//插入操作
				mks[L] = HO.Eval(MK[L], DM, N);			//记录该解的评价值
			}
		}

		//解的更新
		vector<double> bS = HO.bestS(mks, MarkovL);		//寻找当前扰动中的最优值
		double bestEval = bS[0];
		int bestIndex = int(bS[1]);
		NS = MK[bestIndex];
		ns = bestEval;

		T = T * Alpha;
	}

	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << "花费了"
		<< double(duration.count()) * microseconds::period::num / microseconds::period::den
		<< "秒" << endl;

	//输出结果
	double gs0 = 15377.711;
	cout << "最优解为" << bs << endl;
	double e = (bs - gs0) / gs0;
	cout << "误差为" << e * 100.0 << '%' << endl;
	cout << "最优路径为" << endl;
	for (int i = 0; i < N; i++) {
		cout << BS[i] + 1 << '\t';
	}
}