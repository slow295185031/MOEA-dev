function Choose = F_choose(FunctionValue,K)

    [N,M] = size(FunctionValue);
    
    %选出边界点
    Choose = false(1,N);
    [~,Extreme] = min(FunctionValue,[],1);
    Choose(Extreme) = true;
    [~,Extreme] = max(FunctionValue,[],1);
    Choose(Extreme) = true;
    
    %选出剩下的点
    if sum(Choose) > K
        Choosed = find(Choose);
        k = randperm(sum(Choose));
        Choose(Choosed(k(1:sum(Choose)-K))) = false;
    elseif sum(Choose) < K
        Distance = zeros(N)+inf;
        for i = 1 : N-1
            for j = i+1 : N
                Distance(i,j) = norm(FunctionValue(i,:)-FunctionValue(j,:));
                Distance(j,i) = Distance(i,j);
            end
        end
        while sum(Choose) < K
            Remain = find(~Choose);
            [~,x] = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
        end
    end
end

