function [DistanceDirect,DistanceMaped] = F_distance(Fun1,Fun2)
%计算两两个体向量间的距离(越大越好)

    [N1,M] = size(Fun1);
    N2 = size(Fun2,1);
    MFun1 = Fun1./repmat(sum(Fun1.^1,2).^1,1,M);
    MFun2 = Fun2./repmat(sum(Fun2.^1,2).^1,1,M);
    
    DistanceDirect = zeros(N1,N2);
    DistanceMaped = zeros(N1,N2);
    for i = 1 : N1
        DistanceDirect(i,:) = (sum((repmat(Fun1(i,:),N2,1)-Fun2).^2,2).^0.5)';
        DistanceMaped(i,:) = (sum((repmat(MFun1(i,:),N2,1)-MFun2).^2,2).^0.5)';
    end
end

