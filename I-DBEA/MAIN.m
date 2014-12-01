%I-DBEA
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,p1,p2] = P_settings('I-DBEA',Problem,M);
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化向量
    Evaluations = Generations*N;
    [N,W] = F_weight(p1,p2,M);
    W(W==0) = 0.000001;
    for i = 1 : N
        W(i,:) = W(i,:)./norm(W(i,:));
    end
    Generations = floor(Evaluations/N);
    
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    z = min(FunctionValue);
    a = F_intercept(FunctionValue);
    
    %开始迭代
    for Gene = 1 : Generations
        %对每个个体执行操作
        for i = 1 : N
            %产生子代
            Offspring = P_generator([Population(i,:);Population(randi([1,N]),:)],Boundary,Coding,1);
            OffFunValue = P_objective('value',Problem,M,Offspring);
            
            %判断子代是否被支配
            if any(sum(FunctionValue<=repmat(OffFunValue,N,1),2)==M)
                continue;
            end
            
            %更新个体
            for j = randperm(N)
                ScaledFun = (FunctionValue(j,:)-z)./(a-z);
                ScaledOffFun = (OffFunValue-z)./(a-z);
                d1_old = sum(W(j,:).*ScaledFun);
                d2_old = norm(ScaledFun-d1_old*W(j,:));
                d1_new = sum(W(j,:).*ScaledOffFun);
                d2_new = norm(ScaledOffFun-d1_new*W(j,:));
                if d2_new < d2_old || d2_new == d2_old && d1_new < d1_old
                    %更新当前个体
                    Population(j,:) = Offspring;
                    FunctionValue(j,:) = OffFunValue;
                    a = F_intercept(FunctionValue);
                    break;
                end
            end
            
            %更新最优理想点
            z = min(z,OffFunValue);
        end

        clc;fprintf('I-DBEA,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%----------------------------------------------------------------------------------------- 
%生成结果
    P_output(Population,toc,'I-DBEA',Problem,M,Run);
end