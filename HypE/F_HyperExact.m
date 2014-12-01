function f = F_HyperExact(points,bounds,k)
%采用精确的递归方法来计算个体的HV贡献度

    Ps = size(points,1); 
    if k < 0
        k = Ps; 
    end
    actDim = size(points,2);
    if length(bounds) == 1
        bounds = repmat(bounds,actDim,1);
    end
    pvec = 1:size(points,1);
    alpha = zeros(1,k);
    for i = 1 : k 
        j = 1:i-1; 
        alpha(i) = prod((k-j)./(Ps-j))./i;
    end
    f = hypesub(size(points,1),points,actDim,bounds,pvec,alpha,k); 
end

function h = hypesub(l,A,actDim,bounds,pvec,alpha,k) 
    h = zeros(1,l); 
    [S,i] = sortrows(A,actDim); 
    pvec = pvec(i); 
    for i = 1 : size(S,1) 
        if i < size(S,1) 
            extrusion = S(i+1,actDim)-S(i,actDim); 
        else
            extrusion = bounds(actDim)-S(i,actDim);
        end
        if actDim == 1
            if i > k
                break; 
            end
            if alpha >= 0
                h(pvec(1:i)) = h(pvec(1:i))+extrusion*alpha(i); 
            end
        elseif extrusion > 0
            h = h+extrusion*hypesub(l,S(1:i,:),actDim-1,bounds,pvec(1:i),alpha,k); 
        end
    end
end

