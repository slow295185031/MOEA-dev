function [Offspring,t] = F_generator(p1,p2,MaxValue,MinValue,t)
%交叉,变异产生一个子代

    D = length(p1);
    
    %遗传操作参数
    ProM = 1/D;     %变异概率
    DisC = 20;     	%交叉参数
    DisM = 20;     	%变异参数
    tic
    %模拟二进制交叉
    beta = zeros(1,D);
    miu = rand(1,D);
    beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
    beta(miu>0.5) = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
    beta = beta.*(-1).^randi([0,1],1,D);
    Offspring = (p1+p2)./2+beta.*(p1-p2)./2;
    t(1)=t(1)+toc;
    tic
    %多项式变异
    k = rand(1,D);
    miu = rand(1,D);
    Temp = (k<=ProM & miu<0.5);
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
    Temp = (k<=ProM & miu>=0.5);
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));        
    t(2)=t(2)+toc;
    %越界处理
    Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
    Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
end