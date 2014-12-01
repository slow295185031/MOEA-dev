function [FrontValue,MaxFront] = P_sort(FunctionValue,Operation)
% 进行非支配排序
% 输入: FunctionValue, 待排序的种群(目标空间)
%       Operation,     可指定仅排序第一个面,排序前一半个体,或是排序所有的个体, 默认为排序所有的个体
% 输出: FrontValue, 排序后的每个个体所在的前沿面编号, 未排序的个体前沿面编号为inf
%       MaxFront,   排序的最大前面编号

    if nargin < 2
        Kind = 1;
    elseif strcmp(Operation,'half')
        Kind = 2;
    elseif strcmp(Operation,'first')
        Kind = 3;
    else
        Kind = 1;
    end
	[N,M] = size(FunctionValue);
    
    MaxFront = 0;
    cz = zeros(1,N);
    FrontValue = zeros(1,N)+inf;
    [FunctionValue,Rank] = sortrows(FunctionValue);
    while (Kind==1 && sum(cz)<N) || (Kind==2 && sum(cz)<N/2) || (Kind==3 && MaxFront<1)
        MaxFront = MaxFront+1;
        d = cz;
        for i = 1 : N
            if ~d(i)
                for j = i+1 : N
                    if ~d(j)
                        k = 1;
                        for m = 2 : M
                            if FunctionValue(i,m) > FunctionValue(j,m)
                                k = 0;
                                break;
                            end
                        end
                        if k == 1
                            d(j) = 1;
                        end
                    end
                end
                FrontValue(Rank(i)) = MaxFront;
                cz(i) = 1;
            end
        end
    end
end


