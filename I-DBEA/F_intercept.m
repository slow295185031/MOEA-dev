function a = F_intercept(FunctionValue)
%计算各维的截距

    [N,M] = size(FunctionValue);

    %找出边界点
    [~,Choosed(1:M)] = min(FunctionValue,[],1);
    L2NormABO = zeros(N,M);
    for i = 1 : M
    	L2NormABO(:,i) = sum(FunctionValue(:,[1:i-1,i+1:M]).^2,2);
    end
    [~,Choosed(M+1:2*M)] = min(L2NormABO,[],1);
    [~,Extreme] = max(FunctionValue(Choosed,:),[],1);
    Extreme = unique(Choosed(Extreme));
    
    %计算截距
    if length(Extreme) < M
        a = max(FunctionValue,[],1);
    else
        Hyperplane = FunctionValue(Extreme,:)\ones(M,1);
        a = 1./Hyperplane';
    end
end

