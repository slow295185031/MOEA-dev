function Population = F_generator(Population,Boundary,FunctionValue,miu,delta,Problem,M)
%个体局部搜索并生成新的种群

    [N,D] = size(Population);
    MaxValue = Boundary(1,:);
    MinValue = Boundary(2,:);
    
    %对每个个体局部搜索出2D个子代
    for i = 1 : N
        for d = 1 : D
            c = miu+delta*randn;
            k = randperm(N);
            w = repmat(Population(i,:),2,1);
            w(1,d) = w(1,d)+c*(Population(k(1),d)-Population(k(2),d));
            w(2,d) = w(2,d)-c*(Population(k(1),d)-Population(k(2),d));
            for j = 1 : 2
                w(j,d) = max(min(w(j,d),MaxValue(d)),MinValue(d));
            end
            wFunValue = P_objective('value',Problem,M,w);
            
           	for j = 1 : 2
                k(j) = any(wFunValue(j,:)<FunctionValue(i,:))-any(wFunValue(j,:)>FunctionValue(i,:));
            end
            if k(1) == -1 && k(2) == -1
                continue;
            elseif k(1) > k(2)
                k = 1;
            elseif k(1) < k(2)
                k = 2;
            else
                k = randi([1,2]);
            end
            Population(i,d) = w(k,d);
            FunctionValue(i,:) = wFunValue(k,:);
        end
    end
end

