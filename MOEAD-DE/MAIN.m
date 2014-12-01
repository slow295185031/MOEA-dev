%MOEA/D-DE
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,H,delta,nr] = P_settings('MOEAD-DE',Problem,M);
    A = 1;	%1.采用切比雪夫方法 2.采用PBI方法
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化向量
    Evaluations = Generations*N;
    [N,W] = F_weight(H,M);
    W(W==0) = 0.000001;
    T = floor(N/10);
    Generations = floor(Evaluations/N);

    %邻居判断
    B = zeros(N);
    for i = 1 : N-1
        for j = i+1 : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %初始化种群
    [Population,Boundary] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    Z = min(FunctionValue);

    %开始迭代
    for Gene = 1 : Generations
        %对每个个体执行操作
        for i = 1 : N
            %选出父母
            if rand < delta
                P = B(i,:);
            else
                P = 1:N;
            end
            k = randperm(length(P));
            
            %产生子代
            Offspring = F_generator(Population(i,:),Population(P(k(1)),:),Population(P(k(2)),:),Boundary);
            OffFunValue = P_objective('value',Problem,M,Offspring);

            %更新最优理想点
            Z = min(Z,OffFunValue);
            
            %更新P中的个体
            c = 0;
            for j = randperm(length(P))
                if c >= nr
                    break;
                end
                if A == 1
                    g_old = max(abs(FunctionValue(P(j),:)-Z).*W(P(j),:));
                    g_new = max(abs(OffFunValue-Z).*W(P(j),:));
                elseif A == 2
                    d1 = abs(sum((FunctionValue(P(j),:)-Z).*W(P(j),:)))/norm(W(P(j),:));
                    g_old = d1+5*norm(FunctionValue(P(j),:)-(Z+d1*W(P(j),:)/norm(W(P(j),:))));               
                    d1 = abs(sum((OffFunValue-Z).*W(P(j),:)))/norm(W(P(j),:));
                    g_new = d1+5*norm(OffFunValue-(Z+d1*W(P(j),:)/norm(W(P(j),:))));
                end               
                if g_new < g_old
                    %更新当前向量的个体
                    Population(P(j),:) = Offspring;
                    FunctionValue(P(j),:) = OffFunValue;
                    c = c+1;
                end
            end
        end
        
        clc;fprintf('MOEA/D-DE,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%----------------------------------------------------------------------------------------- 
%生成结果
    P_output(Population,toc,'MOEAD-DE',Problem,M,Run);
end