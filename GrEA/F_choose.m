function Choose = F_choose(FunctionValue,K,div)
%环境选择,利用格子选出指定个数的最佳个体

    [N,M] = size(FunctionValue);

    %计算格子坐标
    fmax = max(FunctionValue);
    fmin = min(FunctionValue);
    lb = fmin-(fmax-fmin)/2/div;
    ub = fmax+(fmax-fmin)/2/div;
    d = (ub-lb)/div;
    lb = repmat(lb,N,1);
    d = repmat(d,N,1);
    GLoc = floor((FunctionValue-lb)./d);
    
    %计算GR,GCD,GCPD,GD
    GR = sum(GLoc,2);
    GCD = zeros(1,N);
    GCPD = sqrt(sum(((FunctionValue-(lb+GLoc.*d))./d).^2,2));
    GD = zeros(N)+inf;
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));
            GD(j,i) = GD(i,j);
        end
    end

    %判断两两个体的格子支配关系
    G = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(GLoc(i,:)<GLoc(j,:))-any(GLoc(i,:)>GLoc(j,:));
            if k == 1
                G(i,j) = true;
            elseif k == -1
                G(j,i) = true;
            end
        end
    end

    %环境选择
    Choose = true(1,N);
    while sum(Choose) > N-K      
        %选出剩下的最好的一个
        CanBeChoose = find(Choose);
        temp = find(GR(CanBeChoose)==min(GR(CanBeChoose)));
        temp2 = find(GCD(CanBeChoose(temp))==min(GCD(CanBeChoose(temp))));
        [~,q] = min(GCPD(CanBeChoose(temp(temp2))));
        q = CanBeChoose(temp(temp2(q)));
        Choose(q) = false;      
        %更新GCD
        GCD = GCD+max(M-GD(q,:),0);
        %更新GR
        Eq = GD(q,:)==0 & Choose;
        Gq = G(q,:) & Choose;
        NGq = Choose.*(1-Gq);
        Nq = GD(q,:)<M & Choose;
        GR(Eq) = GR(Eq)+M+2;
        GR(Gq) = GR(Gq)+M;
        PD = zeros(N,1);
        for p = find((Nq.*NGq).*(1-Eq))
            if PD(p) < M-GD(q,p)
                PD(p) = M-GD(q,p);
                Gp = G(p,:) & Choose;
                for r = find(Gp.*(1-(Gq+Eq)))
                    if PD(r) < PD(p)
                        PD(r) = PD(p);
                    end
                end
            end
        end
        pp = logical(NGq.*(1-Eq));       
        GR(pp) = GR(pp)+PD(pp);
    end
    Choose =~ Choose;
end

