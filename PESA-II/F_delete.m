function Del = F_delete(FunctionValue,K,div)
    
    [N,M] = size(FunctionValue);

    %计算格子坐标
    fmax = max(FunctionValue);
    fmin = min(FunctionValue);
    d = (fmax-fmin)/div;
    fmin = repmat(fmin,N,1);
    d = repmat(d,N,1);
    GLoc = floor((FunctionValue-fmin)./d);
    GLoc(GLoc>=div) = div-1;

    %确定每个个体所在的格子
    UniqueGLoc = sortrows(unique(GLoc,'rows'));
    [~,Site] = ismember(GLoc,UniqueGLoc,'rows');

    %计算每个格子的拥挤度
    Temp = sortrows(tabulate(Site));
    CrowdG = Temp(:,2);

    %每次删去一个
    Del = false(1,N);
    while sum(Del) < K
        %选择最挤的格子(若同时有多个,则随机选择一个)
        maxGrid = find(CrowdG==max(CrowdG));
        Temp = randi([1,length(maxGrid)]);
        Grid = maxGrid(Temp);
        %从中随机删去一个个体
        InGrid = find(Site==Grid);
        Temp = randi([1,length(InGrid)]);
        p = InGrid(Temp);
        Del(p) = true;
        CrowdG(Grid) = CrowdG(Grid)-1;
        Site(p) = NaN;
    end
end

