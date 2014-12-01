function Start(Algorithm,Problem,Objectives,Run)
%批量运行

    if nargin < 4;Run = 1;end
    if ~iscell(Algorithm);Algorithm = {Algorithm};end
    if ~iscell(Problem);Problem = {Problem};end
    for R = Run
        for A = Algorithm
            a = cell2mat(A);
            if exist(a,'dir') == 7
                addpath(a);
            else
                error(['算法 ',a,' 不存在.']);
            end
            for P = Problem
                for M = Objectives
                    MAIN(cell2mat(P),M,R);
                end
            end
        end
    end
    sp=actxserver('SAPI.SpVoice');
    sp.Speak('job finished');
end

