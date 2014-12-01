function [Archive,ArchiveFunValue] = F_archive(C,CFunValue,Archive,ArchiveFunValue,epsilon)
%尝试将个体C加入外部种群

    [NoA,M] = size(ArchiveFunValue);

    %计算格子位置
    CGrid = floor((CFunValue-min(ArchiveFunValue,[],1))/epsilon);
    AGrid = floor((ArchiveFunValue-repmat(min(ArchiveFunValue,[],1),NoA,1))/epsilon);
    
    %判断ε支配关系
    Dominated = find(sum(AGrid-repmat(CGrid,NoA,1)<=0,2)==M,1);
    if isempty(Dominated)
        Dominate = find(sum(repmat(CGrid,NoA,1)-AGrid<=0,2)==M);
        if ~isempty(Dominate)
            Archive(Dominate,:) = [];
            ArchiveFunValue(Dominate,:) = [];
            Archive = [Archive;C];
            ArchiveFunValue = [ArchiveFunValue;CFunValue];
        else
            SameGrid = find(ismember(AGrid,CGrid,'rows'),1);
            if isempty(SameGrid)
                Archive = [Archive;C];
                ArchiveFunValue = [ArchiveFunValue;CFunValue];
            else
                B = CGrid*epsilon+min(ArchiveFunValue,[],1);
                if norm(CFunValue-B) < norm(ArchiveFunValue(SameGrid,:)-B)
                    Archive(SameGrid,:) = C;
                    ArchiveFunValue(SameGrid,:) = CFunValue;
                end
            end
        end
    end
end

