function MatingPool = F_mating(Population,FunctionValue,div)
%交配池选择

    [N,D] = size(Population);

    %计算格子坐标
    fmax = max(FunctionValue);
    fmin = min(FunctionValue);
    lb = fmin-(fmax-fmin)/2/div;
    ub = fmax+(fmax-fmin)/2/div;
    d = (ub-lb)/div;
    lb = repmat(lb,N,1);
    d = repmat(d,N,1);
    GLoc = floor((FunctionValue-lb)./d);
    
    %计算GD
    GD = zeros(N)+inf;
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));
            GD(j,i) = GD(i,j);
        end
    end
    
    %计算GCD
    GD = max(size(FunctionValue,2)-GD,0);
    GCD = sum(GD,2);
    
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
            Domi = any(FunctionValue(p,:)<FunctionValue(q,:))-any(FunctionValue(p,:)>FunctionValue(q,:));
            GDomi = any(GLoc(p,:)<GLoc(q,:))-any(GLoc(p,:)>GLoc(q,:));
            if Domi == 1 || GDomi == 1
                k(j) = p;
            elseif Domi == -1 || GDomi == -1
                k(j) = q;
            elseif GCD(p) < GCD(q)
                k(j) = p;
            elseif GCD(p) > GCD(q)
                k(j) = q;
            elseif rand < 0.5
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

