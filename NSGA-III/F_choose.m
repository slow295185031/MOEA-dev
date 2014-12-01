function Choose = F_choose(FunctionValue1,FunctionValue2,K,Z)
%环境选择

    FunctionValue = [FunctionValue1;FunctionValue2];
    [N,M] = size(FunctionValue);
    N1 = size(FunctionValue1,1);
    N2 = size(FunctionValue2,1);
    NoZ = size(Z,1);

    %目标函数值归一化
    Zmin = min(FunctionValue,[],1);	%每维最小值
    Extreme = zeros(1,M);           %每维的边界点
    w = zeros(M)+0.000001+eye(M);
    for i = 1 : M                   %找出每维边界点
        [~,Extreme(i)] = min(max(FunctionValue./repmat(w(i,:),N,1),[],2));
    end
    Hyperplane = FunctionValue(Extreme,:)\ones(M,1);	%计算超平面
    a = 1./Hyperplane;             	%计算每维的截距
    if any(isnan(a))
        a = max(FunctionValue,[],1)';
    end
    FunctionValue = (FunctionValue-repmat(Zmin,N,1))./(repmat(a',N,1)-repmat(Zmin,N,1));	%归一化
    
    %将每个个体关联到某个参考点
    Distance = zeros(N,NoZ);        %计算每个个体到每个参考点向量的距离
    normZ = sum(Z.^2,2).^0.5;
    normF = sum(FunctionValue.^2,2).^0.5;
    for i = 1 : N
        normFZ = sum((repmat(FunctionValue(i,:),NoZ,1)-Z).^2,2).^0.5;
        for j = 1 : NoZ
            S1 = normF(i);
            S2 = normZ(j);
            S3 = normFZ(j);
            p = (S1+S2+S3)/2;
            Distance(i,j) = 2*sqrt(p*(p-S1)*(p-S2)*(p-S3))/S2;
        end
    end
    [d,pi] = min(Distance',[],1);   %离每个个体最近的参考点向量
    
    %计算每个参考点所关联的除最后一个面的个体数
    rho = zeros(1,NoZ);
    for i = 1 : N1
        rho(pi(i)) = rho(pi(i))+1;
    end
    
    %环境选择
    Choose = false(1,N2);       %最终选出的个体编号
    Zchoose = true(1,NoZ);      %剩余未删除的参考点
    k = 1;
    while k <= K
        Temp = find(Zchoose);
        [~,j] = min(rho(Temp));
        j = Temp(j);
        I = find(Choose==0 & pi(N1+1:end)==j);
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi([1,length(I)]);
            end
            Choose(I(s)) = true;
            rho(j) = rho(j)+1;
            k = k+1;
        else
            Zchoose(j) = false;
        end
    end
end

