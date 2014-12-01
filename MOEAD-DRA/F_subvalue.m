function SubValue = F_subvalue(FunctionValue,W,Z,A)
%计算各向量上的单目标值

    N = size(FunctionValue,1);
    
    SubValue = zeros(1,N);   
    if A == 1
        for i = 1 : N
        	SubValue(i) = max(abs(FunctionValue(i,:)-Z).*W(i,:));
        end
    elseif A == 2
        for i = 1 : N
            d1 = abs(sum((FunctionValue(i,:)-Z).*W(i,:)))/norm(W(i,:));
            SubValue(i) = d1+5*norm(FunctionValue(i,:)-(Z+d1*W(i,:)/norm(W(i,:))));
        end
    end   
end

