function MatingPool = F_mating(Population,FunctionValue,RefPoint,NoSample)
%交配池选择

    [N,D] = size(Population);
    M = size(FunctionValue,2);
    
    %计算个体适应度值
    %3维及以上时采用估计算法
    if M > 2
        k = N;
        F = zeros(1,N);
        alpha = zeros(1,N); 
        for i = 1 : k 
            alpha(i) = prod((k-[1:i-1])./(N-[1:i-1]))./i; 
        end
        Fmin = min(FunctionValue,[],1);
        S = rand(NoSample,M).*repmat(RefPoint-Fmin,NoSample,1)+repmat(Fmin,NoSample,1);
        PdS = false(N,NoSample);
        dS = zeros(1,NoSample);
        for i = 1 : N
            x = sum(repmat(FunctionValue(i,:),NoSample,1)-S<=0,2)==M;
            PdS(i,x) = true;
            dS(x) = dS(x)+1;
        end
        for i = 1 : N
            F(i) = sum(alpha(dS(PdS(i,:))));
        end  
    %2维时采用精确算法
    else
        F = F_HyperExact(FunctionValue,RefPoint,3);
    end
    
    %二元联赛选择
    MatingPool = zeros(N,D);
    Rank = randperm(N);
    Pointer = 1;
    for i = 1 : 2 : N
        %选择父母
        k = zeros(1,2);
        for j = 1 : 2
            if Pointer >= N
                Rank = randperm(N);
                Pointer = 1;
            end
            p = Rank(Pointer);
            q = Rank(Pointer+1);
            if F(p) > F(q)
                k(j) = p;
            else
                k(j) = q;
            end
            Pointer = Pointer+2;
        end
        MatingPool(i,:) = Population(k(1),:);
        MatingPool(i+1,:) = Population(k(2),:);
    end
end

