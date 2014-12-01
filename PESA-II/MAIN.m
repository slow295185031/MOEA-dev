%PESA-II
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,div] = P_settings('PESA-II',Problem,M);
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    
    %开始迭代
    for Gene = 1 : Generations
        %产生子代
        MatingPool = F_mating(Population,FunctionValue,N,div);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = P_objective('value',Problem,M,Population);
        FrontValue = P_sort(FunctionValue,'first');
        
        %选出非支配个体
        Next = FrontValue==1;
        
        %裁剪
        if sum(Next) > N
            Del = F_delete(FunctionValue(Next,:),sum(Next)-N,div);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        
        %下一代种群
        Population = Population(Next,:);
        FunctionValue = FunctionValue(Next,:);
        
        clc;fprintf('PESA-II,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%-----------------------------------------------------------------------------------------     
%生成结果
    P_output(Population,toc,'PESA-II',Problem,M,Run);
end
