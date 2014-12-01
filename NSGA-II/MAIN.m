%NSGA-II
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N] = P_settings('NSGA-II',Problem,M);
	FrontColor = ['ro';'bo';'go';'ko';'r^';'g^';'b^';'k^';'r*';'g*';'b*';'k*';'r>';'g>';'b>';'k>';'r<';'g<';'b<';'k<';'rv';'gv';'bv';'kv'];
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    %This is for test -----Start
    Seq = zeros(size(FunctionValue));
    for i = 1 : M
        [~,tmp] = sort(FunctionValue(:,i));
        [~,Seq(:,i)] = sort(tmp);
    end
    %this is for test -----Start
    FrontValue = P_sort(FunctionValue);
    CrowdDistance = F_distance(FunctionValue,FrontValue);
    
    %开始迭代
    for Gene = 1 : Generations    
        %产生子代
        MatingPool = F_mating(Population,FrontValue,CrowdDistance);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = P_objective('value',Problem,M,Population);
        [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
        CrowdDistance = F_distance(FunctionValue,FrontValue);

        
        %选出非支配的个体        
        Next = zeros(1,N);
        NoN = numel(FrontValue,FrontValue<MaxFront);
        Next(1:NoN) = find(FrontValue<MaxFront);
        
        %选出最后一个面的个体
        Last = find(FrontValue==MaxFront);
        [~,Rank] = sort(CrowdDistance(Last),'descend');
        Next(NoN+1:N) = Last(Rank(1:N-NoN));
        
        %下一代种群
        Population = Population(Next,:);
        FrontValue = FrontValue(Next);
        CrowdDistance = CrowdDistance(Next);
        
		FunctionValue = P_objective('value',Problem,M,Population);
		cla;
		for i = 1 : MaxFront
			FrontCurrent = find(FrontValue==i);
			P_draw(FunctionValue(FrontCurrent,:),FrontColor(i,:));
			pause();
		end
	%	pause();
        %This is for test
%        if Gene == floor(Generations / 2)
%            %This is for test -----Start
%            Seq2 = zeros(size(FunctionValue));
%            for i = 1 : M
%                [~,tmp] = sort(FunctionValue(:,i));
%                [~,Seq2(:,i) ]= sort(tmp);
%            end
%            %this is for test -----Start
%        end
        %End Test
        clc;fprintf('NSGA-II,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
	%-----------------------------------------------------------------------------------------     
    %Seq3 = zeros(size(FunctionValue));
    %for i = 1 : M
    %    [~,tmp] = sort(FunctionValue(:,i));
    %    [~,Seq3(:,i)] = sort(tmp);
    %end
%生成结果
    P_output(Population,toc,'NSGA-II',Problem,M,Run);
end
