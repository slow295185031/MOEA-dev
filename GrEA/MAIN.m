%GrEA
function MAIN(Problem,M,Run)
clc;format compact;tic;
 %-----------------------------------------------------------------------------------------
 %参数设定
    [Generations,N,div] = P_settings('GrEA',Problem,M);
 %-----------------------------------------------------------------------------------------
 %算法开始
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);

    %开始迭代
    for Gene = 1 : Generations
        %产生子代
        MatingPool = F_mating(Population,FunctionValue,div);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = P_objective('value',Problem,M,Population);
        [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
        
        %选出前若干面的个体     
        Next = zeros(1,N);
        NoN = numel(FrontValue,FrontValue<MaxFront);
        Next(1:NoN) = find(FrontValue<MaxFront);
        
        %选出最后一个面的个体
        Last = find(FrontValue==MaxFront);
        Choose = F_choose(FunctionValue(Last,:),N-NoN,div);
        Next(NoN+1:N) = Last(Choose);
        
        %下一代种群
        Population = Population(Next,:);
        FunctionValue = FunctionValue(Next,:);
        
        cla;
        P_draw(FunctionValue);
        pause();
        
        clc;fprintf('GrEA,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%-----------------------------------------------------------------------------------------     
%生成结果
    P_output(Population,toc,'GrEA',Problem,M,Run);
end