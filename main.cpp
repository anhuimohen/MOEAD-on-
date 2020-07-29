#include"question.h"
#include"evaluation.h"
#include"moead.h"
#include<ctime>
vector< vector<double> > ws;
#define DD 200															//��������
int main()
{
	srand((unsigned int)time(NULL));
	///////////////////////////////////
	//  ���ľ��߿ռ��Сʱ �Ƚ��г�ʼ��д���ļ�
	/*
	decision_space spa;
	spa.init();															//���߿ռ�ĳ�ʼ��
	spa.print();														//���߿ռ��ӡ
	spa.write_true_value();												//�����߿ռ�д��txt�ļ�
	spa.write_normalization_value();
	*/
	//////////////////////////////////
	//  ��ȡ���߿ռ�����
	
	///////////////////////////////////
	//  �㷨IGDֵ�ļ���
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
	pop.init(space);													//��Ⱥ��ʼ��													
	cal_z_min(pop);														//����Z*
	z_min_print();
	cout<<pop.mean_vector(objective, 5);								//���ɾ��ȷֲ�������
	pop.vector_print();
	pop.cal_nei_index();												//�����������ھ���������
	pop.cal_value_g();
	pop.print();
	save_pop_data(pop,1);
	pop.update(space,1);
	int dd = 1;
	do
	{
		cout << "��" << dd << "�ε�����+++++++++++++++++++";
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