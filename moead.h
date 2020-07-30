#pragma once
#include"question.h"
#include <algorithm>
#include <fstream>
#include<math.h>
#include<string>
#define popsize 126                                                //种群规模
#define num_T  10                                                 //向量邻居的个数
#define p_variate 0.1											 //变异概率
#define SBX_n 2													 //分布指标 n越大表示子代和父代更接近
#define apap_x 400												 //自适应变异参数
#define DE_F 0.5												 //DE算子的缩放因子

class individual;
class populatioon;
double z_min[objective];										 //Z* 最小值
vector<int> num_update ;										 //种群中个体的更新次数
vector<double> mean_objective[objective];						 //记录随着进化代数的增加非支配种群目标值的均值
vector<individual> EP;											 //外部种群
double p_select_concrete[dimension][concrete_service_num];		 //记录当前EP中每维每个具体服务出现的概率

class individual
{
public:
	int select[dimension];										 //每一维选择的服务编号
	double value[objective];                                     //个体解的目标值
	void random_select();										 //随机生成解
	void cal_value(decision_space a);
	void print();
};

void individual::random_select()								 //随机生成解
{
	for (int i = 0; i < dimension; i++)
	{
		select[i] = rand() % concrete_service_num;
	}
}

void individual::cal_value(decision_space a)						     //个体初始化
{
	for (int i = 0; i < objective; i++)
	{
		if (ob_cal_mode[i] == 0)								 //0-相加 初始为0
			value[i] = 0;
		if (ob_cal_mode[i] == 1)                                 //1-想乘 初始为1
			value[i] = 1;
	}
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < objective; j++)
		{
			if (ob_cal_mode[j] == 0)
				value[j] += a.service_unit[i][select[i]].attribute_value[j][1];
			if (ob_cal_mode[j] == 1)							
				value[j] = value[j] * a.service_unit[i][select[i]].attribute_value[j][1]; 
		}
	}
}

void individual::print()                                         //打印个体
{
	cout << "个体每维选择的是：";
	for (int i = 0; i < dimension; i++)
	{
		cout << setw(5) << select[i];
	}
	cout << "目标值是：";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(15) << value[i];
	}
}

class population
{
public:
	individual pop[popsize];
	double value_g[popsize];									 //根据切比雪夫公式计算的目标值 g
	vector< vector<double> > ws;								 //权重向量
	int neighbor_index[popsize][num_T];							 //每向量的邻居向量编号
	int mean_vector(int m, double H);						     //生成均匀分布的向量
	void vector_print();
	void init(decision_space a);
	void print();
	void cal_nei_index();
	void cal_value_g();											 //切比雪夫公式计算的目标值 g
	void update(decision_space space,int dd);					 //主操作：种群更新
	void select_parents(int k, individual& z, individual& x);
};

void population::init(decision_space a)
{
	for (int i = 0; i < objective; i++)							//初始z*_min													
	{
		z_min[i] = RAND_MAX;
	}
	for (int i = 0; i < popsize; i++)
	{
		pop[i].random_select();
		pop[i].cal_value(a);
	}

}

void population::print()
{
	cout << "打印种群：*******************" << endl;
	for (int i = 0; i < popsize; i++)
	{
		pop[i].print();
		cout << "   value_g  " << value_g[i];
		cout << endl;
	}
}

void cal_z_min(population a)                                    //计算各个目标的最小值 即z*
{
	for (int i = 0; i < objective; i++)							//初始z*_min													
	{
		z_min[i] = RAND_MAX;
	}
	for (int i = 0; i < popsize; i++)
	{
		for (int j = 0; j < objective; j++)
		{
			if (a.pop[i].value[j] < z_min[j])
				z_min[j] = a.pop[i].value[j];
		}
	}
}

void z_min_update(individual a)									//利用个体a更新Z*
{
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] < z_min[i])
			z_min[i] = a.value[i];
	}
}

void z_min_print()
{
	cout << "打印各个目标的最小值Z*：";
	for (int i = 0; i < objective; i++)
	{
		cout << setw(15) << z_min[i];
	}
	cout << endl;
}

int population::mean_vector(int m, double H)//在m维目标空间的每维中均匀产生H个向量
{
	/*
	int m; // the number of objectives
	double stepsize;
	double H; // H = 1 / stepsize
	cout << "Please input the number of objectives (m): \n";
	cin >> m;
	cout << "Please input the stepsize (1/H): \n";
	cin >> stepsize;
	H = 1 / stepsize;
	cout << "H = " << H << endl;
	*/
	vector<int> sequence;
	for (unsigned int i = 0; i < H; i++) // the number of zero is (H)
	{
		sequence.push_back(0);
	}
	for (unsigned int i = 0; i < (m - 1); i++) // the number of 1 is (H + m - 1 - (m - 1))
	{
		sequence.push_back(1);
	}

	do
	{
		int s = -1;
		vector<double> weight;
		for (unsigned int i = 0; i < (H + m - 1); i++)
		{
			if (sequence[i] == 1)
			{
				double w = i - s;
				w = (w - 1) / H;
				s = i;
				weight.push_back(w);
			}
		}
		double w = H + m - 1 - s;
		w = (w - 1) / H;
		weight.push_back(w);
		ws.push_back(weight);
	} while (next_permutation(sequence.begin(), sequence.end()));
	ofstream outfile("weight.txt");
	for (unsigned int i = 0; i < ws.size(); i++)
	{
		for (unsigned int j = 0; j < ws[i].size(); j++)
		{
			outfile << ws[i][j] << " ";
		}
		outfile << "\n";
	}
	return ws.size();
}

void population::vector_print()
{
	cout << "vector:" << endl;
	for (int i = 0; i < ws.size(); i++)
	{
		for (int j = 0; j < ws[0].size(); j++)
		{
			cout << ws[i][j] << '\t';
		}
		cout << endl;
	}
}

void population::cal_nei_index()                                  //计算向量的邻居向量的索引
{
	for (int i = 0; i < popsize; i++)
	{
		int index[popsize];
		double distance[popsize];
		for (int j = 0; j < popsize; j++)
		{
			index[j] = j;
			double dis=0;
			for (int k = 0; k < objective; k++)
			{
				dis += pow((ws[i][k] - ws[j][k]), 2);
			}
			distance[j] = sqrt(dis);
		}
		for (int j = 0; j < popsize; j++)                       //将距离按照从小到大排列
		{														//同时索引编号也更换位置
			for (int k = 0; k < popsize - 1; k++)
			{
				if (distance[k] > distance[k + 1])
				{
					double di = distance[k];
					distance[k] = distance[k + 1];
					distance[k + 1] = di;
					int in = index[k];
					index[k] = index[k + 1];
					index[k + 1] = in;
				}
			}
		}
		for (int j = 0; j < num_T; j++)
		{
			neighbor_index[i][j] = index[j];
		}
	}
	////////////////////////////////////////////
	// 输出每个向量的邻居向量的索引
	for (int i = 0; i < popsize; i++)
	{
		cout << "第" << i << "个向量的邻居是：";
		for (int j = 0; j < num_T; j++)
		{
			cout << neighbor_index[i][j] << '\t';
		}
		cout << endl;
	}
}

void population::cal_value_g()										//切比雪夫公式计算的目标值 g
{
	for (int i = 0; i < popsize; i++)
	{
		double temporary_g=0;
		for (int j = 0; j < objective; j++)
		{
			double tem = (pop[i].value[j] - z_min[j]) * ws[i][j];
			if (tem > temporary_g)
				temporary_g = tem;
		}
		value_g[i] = temporary_g;
	}
}

void DE(individual& a, individual& b,individual &c)											//DE算子
{
	int x = rand() % dimension;
	for (int i = x; i < dimension; i++)
	{
		int x1 = round(a.select[i] + DE_F * (b.select[i] - c.select[i]));
		int x2 = round(b.select[i] + DE_F * (a.select[i] - c.select[i]));
		///////////////////////////////////
		//  镜像处理  即超出边界取余处理
		if (x1 < 0)
			x1 = 0 + (0 - x1) % (concrete_service_num - 1);
		if (x1 >= concrete_service_num)
			x1 = concrete_service_num - 1 - (x1 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		if (x2 < 0)
			x2 = 0 + (0 - x2) % (concrete_service_num - 1);
		if (x2 >= concrete_service_num)
			x2 = concrete_service_num - 1 - (x2 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		a.select[i] = x2;
		b.select[i] = x1; 
	}
}

void SBX(individual& a, individual& b)											//SBX算子
{
	int x = rand() % dimension;
	double u = rand() % 100 / 100.0;
	double beita;
	if (u <= 0.5)
		beita = pow(2 * u, 1 / (SBX_n + 1));
	else
		beita = pow(1 / (2 - 2 * u), 1 / SBX_n + 1);
	for (int i = x; i < dimension; i++)
	{
		int c1 = round(0.5 * (double(a.select[i] + b.select[i])) - 0.5 * beita * (double(b.select[i] - a.select[i])));
		int c2 = round(0.5 * (double(a.select[i] + b.select[i])) + 0.5 * beita * (double(b.select[i] - a.select[i])));
		///////////////////////////////////
		//  镜像处理  即超出边界取余处理
		if (c1 < 0)
			c1 = 0 + (0 - c1) % (concrete_service_num - 1);
		if (c1 >= concrete_service_num)
			c1 = concrete_service_num - 1 - (c1 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		if (c2 < 0)
			c2 = 0 + (0 - c2) % (concrete_service_num - 1);
		if (c2 >= concrete_service_num)
			c2 = concrete_service_num - 1 - (c2 - (concrete_service_num - 1)) % (concrete_service_num - 1);
		a.select[i] = c2;
		b.select[i] = c1;
	}
}

void adaptive_variation(individual& a, int dd,decision_space space)			//自适应变异  dd进化代数
{
	double p_mu = 0.2 * (1 - exp(-1 * dd / apap_x));
	int flag = 0;															//是否增加突变概率	1-增加
	if (dd > 10)
	{
		int count = 0;
		for (int i = dd-1; i >= dd - 5; i--)
		{
			if (num_update[i] < 5)											//如果连续5代更新次数均<5，就 增加0.1的变异概率
				count++;
		}
		if (count == 5)
			flag = 1;
	}
	if (flag)
		p_mu += 0.1;
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_mu)
			a.select[i] = rand() % concrete_service_num;
	}
	a.cal_value(space);
}

void adaptive_variation2(individual& a, int dd, decision_space space)			//自适应变异  dd进化代数
{
	double p_mu = 0.2 * (1 - exp(-1 * dd / apap_x));
	int flag = 0;															//是否增加突变概率	1-增加
	if (dd > 10)
	{
		int count = 0;
		for (int i = dd - 1; i >= dd - 5; i--)
		{
			if (num_update[i] < 5)											//如果连续5代更新次数均<5，就 增加0.1的变异概率
				count++;
		}
		if (count == 5)
			flag = 1;
	}
	if (flag)
		p_mu += 0.1;
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_mu)
		{
			double p1 = rand() % 100 / 100.0;
			for (int j = 0; j < concrete_service_num; j++)
			{
				if (p1 < p_select_concrete[i][j])
				{
					a.select[i] = j;
					break;
				}
			}
		}
	}
	a.cal_value(space);
}

void individual_cross(individual& a, individual& b)								//个体交叉
{
	int x1 = rand() % dimension;
	int x2;
	do
	{
		x2 = rand() % dimension;
	} while (x1==x2||x1-x2==dimension-1||x2-x1==dimension-1);
	if (x1 > x2)
	{
		int x3 = x1;
		x1 = x2;
		x2 = x3;
	}
	for (int i = x1; i <= x2; i++)
	{
		int x = a.select[i];
		a.select[i] = b.select[i];
		b.select[i] = x;
	}
}

void cal_p_select_concrete()
{
	memset(p_select_concrete, 0, sizeof(p_select_concrete));
	for (int i = 0; i < EP.size(); i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			p_select_concrete[j][EP[i].select[j]]++;
		}
	}
	for (int i = 0; i < dimension; i++)
	{
		int count = 0;
		for (int j = 0; j < concrete_service_num; j++)
		{
			count += p_select_concrete[i][j];
		}
		for (int j = 0; j < concrete_service_num; j++)
		{
			p_select_concrete[i][j] = p_select_concrete[i][j] / (count * 1.0);
		}
		for (int j = 1; j < concrete_service_num; j++)
		{
			p_select_concrete[i][j] += p_select_concrete[i][j - 1];
		}
	}
}

void p_variation(individual& a, decision_space space)						//依据概率进行轮盘赌发变异
{
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_variate)
		{
			double p1 = rand() % 100 / 100.0;
			for (int j = 0; j < concrete_service_num; j++)
			{
				if (p1 < p_select_concrete[i][j])
				{
					a.select[i] = j;
					break;
				}
			}
		}	
	}
	a.cal_value(space);
}

void variation(individual &a,decision_space space)							//个体变异 计算目标值
{
	for (int i = 0; i < dimension; i++)
	{
		if (rand() % 100 / 100.0 < p_variate)
			a.select[i] = rand() % concrete_service_num;
	}
	a.cal_value(space);
}

double cal_y_g(individual a, vector<double> w)								//计算个体在向量w上的目标值value_g
{
	double temporary=0;
	for (int i = 0; i < objective; i++)
	{
		double tem = (a.value[i] - z_min[i]) * w[i];
		if (tem > temporary)
			temporary = tem;
	}
	return temporary;
}

bool operator<(individual a, individual b)						//计算a是否支配b
{
	int falg = 1;												//1-每一个目标值a<=b
	int fa = 0;													//1-存在一个目标值a<b
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] > b.value[i])
		{
			falg = 0;
			break;
		}
	}
	for (int i = 0; i < objective; i++)
	{
		if (a.value[i] < b.value[i])
		{
			fa = 1;
			break;
		}
	}
	if (falg && fa)
		return true;
	return false;
}

void update_ep(individual a)												//更新外部存档种群EP
{
	////////////////////////////////////////////////
	//  首先删除EP中被a支配的解
	for (int i = 0; i < EP.size(); i++)
	{
		if (a < EP[i])
		{
			EP.erase(EP.begin() + i);
			i--;
		}
	}
	//  如果没有解支配a,就将a加入EP
	int flag = 1;															//1-无个体支配a
	for (int i = 0; i < EP.size(); i++)
	{
		int count = 0;														//相同的目标个数
		if (EP[i] < a)
		{
			flag = 0;
			break;
		}
		for (int j = 0; j < objective; j++)
		{
			if (EP[i].value[j] == a.value[j])
				count++;
		}
		if (count == objective)
		{
			flag = 0;
			break;
		}
	}
	if (flag)
		EP.push_back(a);
}

void ep_print()
{
	cout << "打印存档种群EP------------------------";
	for (int i = 0; i < EP.size(); i++)
	{
		cout << "目标值是:";
		for (int j = 0; j < objective; j++)
		{
			cout << setw(15) << EP[i].value[j];
		}
		cout << endl;
	}
}

void save_num_upate()
{
	string str1 = ".\\data\\num_update\\";
	string str2 = "num_update.xls";
	str2 = str1 + str2;
	ofstream outfile(str2);
	for (int i = 0; i < num_update.size(); i++)
	{
		outfile << num_update[i] << endl;
	}
	outfile.close();
}

void cal_mean_objective()													//计算EP种群目标值的均值
{
	double ob[objective] = { 0 };
	for (int i = 0; i < objective; i++)
	{
		for (int j = 0; j < EP.size(); j++)
			ob[i] += EP[j].value[i];
		ob[i] = ob[i] / EP.size();
		mean_objective[i].push_back(ob[i]);
	}
}

void save_mean_objective()													//将得到的目标均值随进化代数的数据保存起来
{
	string str1 = ".\\data\\num_update\\";
	string str2 = "EP_objective_mean_dd.xls";
	str2 = str1 + str2;
	ofstream out(str2);
	for (int i = 0; i < mean_objective[0].size(); i++)
	{
		for (int j = 0; j < objective; j++)
		{
			out << mean_objective[j][i] << '\t';
		}
		out << endl;
	}
	out.close();
}

void population::select_parents(int k,individual &z,individual &x)		
{																			//选择操作，选择个体邻域中较优的个体
	vector<individual> non_dominant;
	for (int i = 0; i < num_T; i++)
	{
		int flag = 1;														//1-第i个邻居个体是非支配的
		for (int j = 0; j < num_T; j++)
		{
			if (i != j)
			{
				if (pop[neighbor_index[k][j]] < pop[neighbor_index[k][i]])
				{
					flag = 0;
					break;
				}
			}
		}
		if (flag)
		{
			int fla = 1;
			for (int q = 0; q < non_dominant.size(); q++)
			{
				int count = 0;
				for (int p = 0; p < objective; p++)
				{
					if (pop[neighbor_index[k][i]].value[p] == non_dominant[q].value[p])
						count++;
				}
				if (count == objective)
				{
					fla = 0; break;
				}
			}
			if (fla)
				non_dominant.push_back(pop[neighbor_index[k][i]]);
		}
	}
	if (non_dominant.size() < 2)
	{
		int x1 = rand() % num_T;											//随机产生个体i的两个邻居编号
		int x2;
		do
		{
			x2 = rand() % num_T;
		} while (x1 == x2 );
		int w = neighbor_index[k][x1];
		int l = neighbor_index[k][x2];
		 z = pop[w];
		 x = pop[l];
	}
	else
	{
		int x1 = rand() % non_dominant.size();
		int x2;
		do
		{
			x2 = rand() % non_dominant.size();
		} while (x1 == x2);
		z = non_dominant[x1];
		x = non_dominant[x2];
	}
}

void population::update(decision_space space,int dd)
{ 
	int number_update = 0;													//每轮个体的更新次数
	for (int i = 0; i < popsize; i++)
	{
		individual a, b, c;
		//////////////////////////////
		//  随机选择邻居中的两个个体
		/*
		int x1 = rand() % num_T;											//随机产生个体i的两个邻居编号
		int x2,x3;
		do
		{
			x2 = rand() % num_T;
			x3 = rand() % num_T;
		} while (x1==x2||x1==x3||x2==x3);
		int k = neighbor_index[i][x1];
		int l = neighbor_index[i][x2];
		int j = neighbor_index[i][x3];
		a = pop[k];
		b = pop[l];
		c = pop[j];
		*/
		//////////////////////////////////
		//  选择邻居中的非支配个体作为父代
		select_parents(i, a, b);
		/*
		double update_rate = 0;												//个体更新频率
		if (dd > 10)
		{
			for (int j = dd-1; j >= dd - 5; j--)
			{
				update_rate += num_update[j];
			}
			update_rate = update_rate / 5;
			if (update_rate < 5)						//如果近10次平均更新频率<5，则采用SBX算子，否则采用DE算子
				SBX(a, b);
			else
				DE(a, b,c);
		}
		else
			DE(a, b,c);
		adaptive_variation(a,dd,space);
		adaptive_variation(b, dd, space);
		*/


		individual_cross(a, b);

		/////////////////////////
		//  依据概率轮盘赌变异  --效果差
		/*
		cal_p_select_concrete();
		p_variation(a, space);
		p_variation(b, space);
		*/
		//////


		variation(a, space);
		variation(b, space);
		z_min_update(a);
		z_min_update(b);
		cal_value_g();
		for (int j = 0; j < num_T; j++)					//对新解在邻居向量的目标值与原个体目标比较，更新个体
		{
			if (cal_y_g(a, ws[neighbor_index[i][j]]) < value_g[neighbor_index[i][j]])
			{
				number_update++;
				pop[neighbor_index[i][j]] = a;
				value_g[neighbor_index[i][j]] = cal_y_g(a, ws[neighbor_index[i][j]]);
			}
			if (cal_y_g(b, ws[neighbor_index[i][j]]) < value_g[neighbor_index[i][j]])
			{
				number_update++;
				pop[neighbor_index[i][j]] = b;
				value_g[neighbor_index[i][j]] = cal_y_g(b, ws[neighbor_index[i][j]]);
			}
		}
		update_ep(a);
		update_ep(b);
	}
	num_update.push_back(number_update);
	cal_mean_objective();
}

void save_pop_data(population a,int dd)
{
	string str1 = ".\\data\\pop\\";
	string str = "_pop_value_.xls";
	string str2=to_string(dd);
	str = str1+str2 + str;
	ofstream outfile(str);
	for (int i = 0; i < popsize; i++)
	{
		for (int j = 0; j < objective; j++)
		{
			outfile << a.pop[i].value[j] << '\t';
		}
		outfile << endl;
	}
	outfile.close();
}

void save_ep_data(int dd)
{
	string str1 = ".\\data\\EP\\";
	string str = "_EP_value_.xls";
	string str2 = to_string(dd);
	str = str1 + str2 + str;
	ofstream outfile(str);
	for (int i = 0; i < EP.size(); i++)
	{
		for (int j = 0; j < objective; j++)
		{
			outfile << EP[i].value[j] << '\t';
		}
		outfile << endl;
	}
	outfile.close();
}