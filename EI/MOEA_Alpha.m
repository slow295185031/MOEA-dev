%测试利用种群编号选择进化
function MOEA_Alpha(Problem,M,Run)
	clc;format compact;tic;
	%初始化参数
	[Generations,N,p1,p2] = P_settings('MOEAD',Problem,M);
	%初始化种群
	[Population,Boundary,Coding] = P_objective('init',Problem,M,N);
	FunctionValue = P_objective('value',Problem,M,Population);
	
	%初始化向量
	Evaluations = Generations * N;
	[N,W] = F_weight(p1,p2,M);
	Generations = floor(Evaluations/N);
	Seq_W = zeros(size(W));
	for i = 1 : M
		[~,tmp] = sort(W(:,i));
		[~,Seq_W(:,i)] = sort(tmp);
	end
	
	%开始迭代
	for Gene = 1 : Generations
		%产生子代
		Offspring = P_generator(Population,Boundary,Coding,N);
		Population = [Population;Offspring];
		FunctionValue = P_objective('value',Problem,M,Population);
		Seq = zeros(size(FunctionValue));
		for i = 1 ： M
			[~,tmp] = sort(FunctionValue(:,i));
			[~,Seq(:,i)] = sort(tmp);
		end
		
	end
end

