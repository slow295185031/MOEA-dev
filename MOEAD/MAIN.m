%MOEA/D
function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------------------
%参数设定
    [Generations,N,p1,p2] = P_settings('MOEAD',Problem,M);
    A = 1;	%1.采用切比雪夫方法 2.采用PBI方法
%-----------------------------------------------------------------------------------------
%算法开始
    %初始化向量
    Evaluations = Generations*N;
    [N,W] = F_weight(p1,p2,M);
    W(W==0) = 0.000001;
    T = floor(N/10);
    Generations = floor(Evaluations/N);

    %邻居判断
    B = zeros(N);
    for i = 1 : N
        for j = i : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %初始化种群
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    Z = min(FunctionValue);

    %开始迭代
    for Gene = 1 : Generations
        %对每个个体执行操作
        for i = 1 : N
            %归一化
            Fmax = max(FunctionValue);
            Fmin = Z;
            FunctionValue = (FunctionValue-repmat(Fmin,N,1))./repmat(Fmax-Fmin,N,1);
            
            %选出父母
            k = randperm(T);
            k = B(i,k(1:2));

            %产生子代
            Offspring = P_generator([Population(k(1),:);Population(k(2),:)],Boundary,Coding,1);
            OffFunValue = P_objective('value',Problem,M,Offspring);
            OffFunValue = (OffFunValue-Fmin)./(Fmax-Fmin);
            
            %更新最优理想点
            Z = min(Z,OffFunValue);

            %更新邻居个体
            for j = 1 : T
                if A == 1
                    g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                elseif A == 2
                    d1 = abs(sum((FunctionValue(B(i,j),:)-Z).*W(B(i,j),:)))/norm(W(B(i,j),:));
                    g_old = d1+5*norm(FunctionValue(B(i,j),:)-(Z+d1*W(B(i,j),:)/norm(W(B(i,j),:))));               
                    d1 = abs(sum((OffFunValue-Z).*W(B(i,j),:)))/norm(W(B(i,j),:));
                    g_new = d1+5*norm(OffFunValue-(Z+d1*W(B(i,j),:)/norm(W(B(i,j),:))));
                end
                if g_new < g_old
                    %更新当前向量的个体
                    Population(B(i,j),:) = Offspring;
                    FunctionValue(B(i,j),:) = OffFunValue;
                end
            end

            %反归一化
            FunctionValue = FunctionValue.*repmat(Fmax-Fmin,N,1)+repmat(Fmin,N,1);

        end
        
        cla;
        P_draw(FunctionValue);
        pause();
        clc;fprintf('MOEA/D,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
    end
%----------------------------------------------------------------------------------------- 
%生成结果
    
    P_output(Population,toc,'MOEAD',Problem,M,Run);
end