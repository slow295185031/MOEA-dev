%MOEA/D-M2M
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,H,S] = P_settings('MOEAD-M2M',Problem,M);
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化向量
    Evaluations = Generations*N;
    [N,W] = F_weight(H,M);
    W(W==0) = 0.000001;
    Generations = floor(Evaluations/N/S);
    
    %初始化种群
    [Population,Boundary] = P_objective('init',Problem,M,N*S);
    FunctionValue = P_objective('value',Problem,M,Population);
    Choose = F_allocation(FunctionValue,W,S);
    Population = Population(Choose,:);

    %开始迭代
    for Gene = 1 : Generations
        %产生子代
        R = zeros(N*S,size(Population,2));
        for n = 1 : N
            R((n-1)*S+1:(n-1)*S+S,:) = F_generator(Population((n-1)*S+1:(n-1)*S+S,:),Boundary);
        end
        
        %更新各子种群
        Q = [R;Population];
        QFunValue = P_objective('value',Problem,M,Q);
        Choose = F_allocation(QFunValue,W,S);
        Population = Q(Choose,:);
        
        clc;fprintf('MOEA/D-M2M,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%----------------------------------------------------------------------------------------- 
%生成结果
    P_output(Population,toc,'MOEAD-M2M',Problem,M,Run);
end