function [Choose,GChoose] = F_choose(FunctionValue,Goal,FrontValue)
%环境选择,先计算适应度值再选出一半的个体与目标点
   
    [N,M] = size(FunctionValue);
    NGoal = size(Goal,1);

    %计算ng
    FdG = false(N,NGoal);
    dG = zeros(1,NGoal);
    for i = 1 : N
        x = sum(repmat(FunctionValue(i,:),NGoal,1)-Goal<=0,2)==M;
        FdG(i,x) = true;
        dG(x) = dG(x)+1;
    end
    
    %计算Fs
    Fs = zeros(1,N);
    for i = 1 : N
        Fs(i) = sum(1./dG(FdG(i,:)));
    end
    
    %计算Fg
    Fg = zeros(1,NGoal);
    for i = 1 : NGoal
        if dG(i) == 0
            Fg(i) = 0.5;
        else
            Fg(i) = 1/(1+(dG(i)-1)/(N-1));
        end
    end   
    
    %选出一半个体
    NF = find(FrontValue==1);
    if length(NF) < N/2
        Fs(NF) = inf;
        [~,Rank] = sort(Fs,'descend');
        Choose = Rank(1:N/2);
    else
        [~,Rank] = sort(Fs(NF),'descend');
        Choose = NF(Rank(1:N/2));
    end
    
    %选出一半目标点
    [~,Rank] = sort(Fg,'descend');
    GChoose = Rank(1:NGoal/2);
end

