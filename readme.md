### MOEAD算法更改目标数时的操作方式

+ 1.需要修改的数值

```cpp
#define objective 5												   //目标数或者属性数

int ob_cal_mode[objective] = { 0,0,1,1,1 };                            //目标值计算方式 0-相加 1-相乘

int attribute_flag[objective] = { 0,0,1,1,1};                     //属性极性，0-消极属性 1-积极属性
```

+ 2.决策空间重新写入
+ 3.均匀向量重新选择大小，更改popsize

### MOEAD改到AWS——MOEAD

首先决策空间不用重新生成，注释掉，直接读取，然后更换种群更新策略中的个体交叉及变异操作

```cpp
/*
		double update_rate = 0;												//个体更新频率
		if (dd > 15)
		{
			for (int j = dd-1; j >= dd - 10; j--)
			{
				update_rate += num_update[j];
			}
			update_rate = update_rate / 10;
			if (update_rate < 5)						//如果近10次平均更新频率<5，则采用DE算子，否则采用SBX算子
				SBX(a, b);
			else
				DE(a, b, c);
		}
		else
			DE(a, b, c);
		adaptive_variation(a,dd,space);
		adaptive_variation(b, dd, space);
		*/
		individual_cross(a, b);
		variation(a, space);
		variation(b, space);
```

变成

```cpp

		double update_rate = 0;												//个体更新频率
		if (dd > 15)
		{
			for (int j = dd-1; j >= dd - 10; j--)
			{
				update_rate += num_update[j];
			}
			update_rate = update_rate / 10;
			if (update_rate < 5)						//如果近10次平均更新频率<5，则采用DE算子，否则采用SBX算子
				SBX(a, b);
			else
				DE(a, b, c);
		}
		else
			DE(a, b, c);
		adaptive_variation(a,dd,space);
		adaptive_variation(b, dd, space);
		
		individual_cross(a, b);
		variation(a, space);
		variation(b, space);
```

###  实验结果对比

```cpp
read_solution(".\\evalution\\aes-8-200_EP_value_.xls",20954,aesmoead_solution);
	read_solution(".\\evalution\\moead-8-200_EP_value_.xls",37318, moead_solution);
	read_solution(".\\evalution\\select-8-200_EP_value_.xls", 49583, select_solution);
	add_solution(aesmoead_solution,non_dominated_solution);
	add_solution(moead_solution, non_dominated_solution);
	add_solution(select_solution, non_dominated_solution);
	cout << cal_IGD(moead_solution, non_dominated_solution) << endl;
	cout << cal_IGD(aesmoead_solution, non_dominated_solution) << endl;
	cout << cal_IGD(select_solution, non_dominated_solution) << endl;

	cout << cal_coverage(moead_solution, select_solution) << endl;
	cout << cal_coverage(select_solution, moead_solution)<<endl;
	cout << cal_coverage(moead_solution, aesmoead_solution) << endl;
	cout << cal_coverage(aesmoead_solution, moead_solution) << endl;
	cout << cal_coverage(aesmoead_solution, select_solution) << endl;
	cout << cal_coverage(select_solution, aesmoead_solution) << endl;
```

输出结果为

0.00352962
0.00748838
0.00217626
0.0909989
0.13945
0.2893
0.0455812
0.0293044
0.319987