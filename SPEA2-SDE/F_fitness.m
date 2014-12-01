function FitnessValue = F_fitness(FunctionValue)

    [N,M] = size(FunctionValue);

    %计算两两个体的支配关系
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(FunctionValue(i,:)<FunctionValue(j,:))-any(FunctionValue(i,:)>FunctionValue(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    %计算S(i)
    S = sum(Dominate,2);
    
    %计算R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end
    
    %计算两两个体的距离
    Distance = zeros(N);
    for i = 1 : N
        ShiftFunctionValue = FunctionValue;
        Temp = repmat(FunctionValue(i,:),N,1);
        Shifted = FunctionValue<Temp;
        ShiftFunctionValue(Shifted) = Temp(Shifted);
        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(FunctionValue(i,:)-ShiftFunctionValue(j,:));
        end
    end
    
    %计算D(i)
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);
    
    %计算适应度值
    FitnessValue = R+D';
end

