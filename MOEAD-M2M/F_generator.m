function NewPopulation = F_generator(Population,Boundary)
%交叉,变异并生成新的子种群

    [N,D] = size(Population);
    
    %遗传操作参数
    ProM = 1/D;     %变异概率
    DisC = 20;     	%交叉参数
    DisM = 20;     	%变异参数
    
    %模拟二进制交叉
    NewPopulation = zeros(N,D);
    for i = 1 : N
        k = randi([1,N-1]);
        if k >= i
            k = k+1;
        end
        beta = zeros(1,D);
        miu = rand(1,D);
        beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
        beta(miu>0.5) = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
        beta = beta.*(-1).^randi([0,1],1,D);
        NewPopulation(i,:) = (Population(i,:)+Population(k,:))./2+beta.*(Population(i,:)-Population(k,:))./2;
    end
    
    %多项式变异
    MaxValue = repmat(Boundary(1,:),N,1);
    MinValue = repmat(Boundary(2,:),N,1);
    k = rand(N,D);
    miu = rand(N,D);
    Temp = (k<=ProM & miu<0.5);
    NewPopulation(Temp) = NewPopulation(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(NewPopulation(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
    Temp = (k<=ProM & miu>=0.5);
    NewPopulation(Temp) = NewPopulation(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-NewPopulation(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));        
    
    %越界处理
    NewPopulation(NewPopulation>MaxValue) = MaxValue(NewPopulation>MaxValue);
    NewPopulation(NewPopulation<MinValue) = MinValue(NewPopulation<MinValue);
end