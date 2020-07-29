#include"question.h"
#include"evaluation.h"
#include"moead.h"
#include<ctime>
vector< vector<double> > ws;
#define DD 200															//迭代次数
int main()
{
	srand((unsigned int)time(NULL));
	///////////////////////////////////
	//  更改决策空间大小时 先进行初始化写入文件
	/*
	decision_space spa;
	spa.init();															//决策空间的初始化
	spa.print();														//决策空间打印
	spa.write_true_value();												//将决策空间写入txt文件
	spa.write_normalization_value();
	*/
	//////////////////////////////////
	//  读取决策空间数据
	
	///////////////////////////////////
	//  算法IGD值的计算
	read_solution(".\\evalution\\moead-5-200_EP_value_02.xls",1976,moead_solution);
	read_solution(".\\evalution\\aes-5-200_EP_value_02.xls",1731, aesmoead_solution);
	add_solution(moead_solution, aesmoead_solution, non_dominated_solution);
	cout << cal_IGD(moead_solution, non_dominated_solution) << endl;
	cout << cal_IGD(aesmoead_solution, non_dominated_solution) << endl;

	cout << cal_coverage(moead_solution, aesmoead_solution) << endl;
	cout << cal_coverage(aesmoead_solution, moead_solution)<<endl;
	


	decision_space space;
	space.read_true_value();
	space.read_normalization_value();
	space.print();
	population pop;
	pop.init(space);													//种群初始化													
	cal_z_min(pop);														//计算Z*
	z_min_print();
	cout<<pop.mean_vector(objective, 5);								//生成均匀分布的向量
	pop.vector_print();
	pop.cal_nei_index();												//计算向量的邻居向量索引
	pop.cal_value_g();
	pop.print();
	save_pop_data(pop,1);
	pop.update(space,1);
	int dd = 1;
	do
	{
		cout << "第" << dd << "次跌代：+++++++++++++++++++";
		//pop.print();
		//ep_print();
		if (dd % 25 == 0)
		{
			save_pop_data(pop, dd);
			save_ep_data(dd);
		}
		pop.update(space,dd);
		dd++;
	} while (dd<DD+1);
	save_num_upate();
	save_mean_objective();
}