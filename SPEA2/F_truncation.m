function Del = F_truncation(FunctionValue,K)
%环境选择

    [N,M] = size(FunctionValue);
    
    %计算两两个体的距离
    Distance = zeros(N)+inf;
    for i = 1 : N-1
        for j = i+1 : N
            Distance(i,j) = norm(FunctionValue(i,:)-FunctionValue(j,:));
            Distance(j,i) = Distance(i,j);
        end
    end
    
    %截断
    Del = false(1,N);
    while sum(Del) < K
        Remain = find(~Del);
        Temp = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end

