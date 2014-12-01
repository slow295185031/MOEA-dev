%PICEA-g
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,NGoal] = P_settings('PICEA-g',Problem,M);
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    Goal = F_Ggene(FunctionValue,NGoal);
    
    %开始迭代
    for Gene = 1 : Generations
        %产生子代
        MatingPool = F_mating(Population);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = P_objective('value',Problem,M,Population);
        Goal = [Goal;F_Ggene(FunctionValue,NGoal)];
        FrontValue = P_sort(FunctionValue,'first');
        
        %环境选择
        [Next,GChoose] = F_choose(FunctionValue,Goal,FrontValue);
        
        %下一代种群
        Population = Population(Next,:);
        Goal = Goal(GChoose,:);

        clc;fprintf('PICEA-g,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%-----------------------------------------------------------------------------------------     
%生成结果
    P_output(Population,toc,'PICEA-g',Problem,M,Run);
end