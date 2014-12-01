function [Generations,N,varargout] = P_settings(Algorithm,Problem,M)
% ����ָ���㷨��ָ�������µĲ�������, ����������,��Ⱥ��С,�㷨�ض������
% ����: Algorithm,  �㷨���
%       Problem,    �����������
%       M,          ��������ά��
% ���: Generations, ������
%       N,           ��Ⱥ��С
%       varargout,   �㷨�ض��Ĳ����; ��ͬ���㷨�Ĳ��������ܲ�ͬ

    [Generations,N] = set_problem(Problem,M);
    Parameter = set_algorithm(Algorithm,Problem,M);
    varargout = num2cell(Parameter);
end

function [Generations,N] = set_problem(Problem,M)
    k = find(~isstrprop(Problem,'digit'),1,'last');
    D = str2double(Problem(k+1:end));
    Problem = Problem(1:k);
    switch Problem
        case 'DTLZ'
            if M < 2 || M > 10
                error('ֻ֧��2~10ά��DTLZ����');
            end
            if D < 1 || D > 7
                error('ֻ֧��DTLZ1~7����');
            end
            Generations = [700 250 1000 250 250 250 250];
            Generations = Generations(D);
        case 'WFG'
            if M < 2 || M > 10
                error('ֻ֧��2~10ά��WFG����');
            end
            if D < 1 || D > 9
                error('ֻ֧��WFG1~9����');
            end
            Generations = [1000 700 250 250 250 250 250 250 250];
            Generations = Generations(D);
        case 'ZDT'
            if M ~= 2
                error('ֻ֧��2ά��ZDT����');
            end
            if D < 1 || D > 6
                error('ֻ֧��ZDT1~6����');
            end
            Generations = 300;
        case 'UF'
            if D >= 1 && D <= 7 && M ~= 2
                error('ֻ֧��2ά��UF1~7����');
            end
            if D >= 8 && D <= 10 && M ~= 3
                error('ֻ֧��3ά��UF8~10����');
            end
            if D < 1 || D > 10
                error('ֻ֧��UF1~10����');
            end
            Generations = 600;
        case 'MOP'
            if M ~= 3
                error('ֻ֧��3ά��MOP����');
            end
            if D < 1 || D > 4
                error('ֻ֧��MOP1~4����');
            end
            Generations = 300;
        case 'MOKP'
            if M < 2 || M > 4
                error('ֻ֧��2~4ά��MOKP����');
            end
            if D ~= 250 && D ~= 500 && D ~= 750
                error('ֻ֧��MOKP250,500,750����');
            end
            Generations = 1000;
        otherwise
            error(['����',Problem,'������.']);
    end
    N = [100 105 120 126 132 112 156 90 275];
    N = N(M-1);
end

function Parameter = set_algorithm(Algorithm,Problem,M)
    Parameter = NaN;
    k = find(~isstrprop(Problem,'digit'),1,'last');
    D = str2double(Problem(k+1:end));
    Problem = Problem(1:k);
    switch Algorithm
        case 'A'
        case 'e-MOEA'
            epsilon = 0.06;
            Parameter(1) = epsilon;
        case 'GrEA'
            switch Problem
                case 'DTLZ'
                    div = [55 16 10 10 10  0 10  0 11
                           45 15 10  9  8  0  7  0  8
                           45 17 11 11 11  0 10  0 11
                           55 15 10  9  8  0  7  0  8
                           55 40 35 29 14  0 11  0 11
                           55 33 36 24 50  0 50  0 50
                           16 15  9  8  6  0  5  0  4];
                    div = div(D,M-1);
                case 'WFG'
                    div = [45 10  8  9  9  0  7  0 10
                           45 12 11 11 11  0 11  0 11
                           55 22 18 18 18  0 16  0 22
                           45 15 10  9  9  0  7  0  8
                           45 15 10  9  9  0  7  0  8
                           45 15 10  9  9  0  7  0  8
                           45 15 10  9  9  0  7  0  8
                           45 15 10  9  9  0  7  0  8
                           45 15 10  9  9  0  7  0  8];
                    div = div(D,M-1);
                otherwise
                    div = [45 15 10  9  9  0  7  0  8];
                    div = div(M-1);
            end
            Parameter(1) = div;
        case 'HypE'
            NoSample = 8888;
            Parameter(1) = NoSample;
        case 'I-DBEA'
            p1 = [99 13  7  5  3  3  3  2  2];
            p2 = [ 0  0  0  0  3  2  1  2  2];
            p1 = p1(M-1);
            p2 = p2(M-1);
            Parameter(1) = p1; Parameter(2) = p2;
        case 'KnEA'
            switch Problem
                case 'DTLZ'
                    rate = [.6 .6 .6 .2 .2  0 .1  0 .1
                            .6 .5 .5 .5 .5  0 .5  0 .5
                            .6 .6 .4 .3 .2  0 .1  0 .1
                            .6 .5 .5 .5 .5  0 .5  0 .5
                            .6 .6 .5 .5 .5  0 .3  0 .3
                            .6 .6 .5 .4 .4  0 .3  0 .3
                            .6 .6 .6 .6 .6  0 .5  0 .6];
                    rate = rate(D,M-1);
                case 'WFG'
                    rate = [.6 .5 .5 .5 .5  0 .5  0 .5];
                    rate = rate(M-1);
                otherwise
                    rate = [.6 .5 .5 .5 .5  0 .5  0 .5];
                    rate = rate(M-1);
            end
            Parameter(1) = rate;
        case 'MOEAD'
            p1 = [99 13  7  5  4  3  3  2  3];
            p2 = [ 0  0  0  0  1  2  2  2  2];
            p1 = p1(M-1);
            p2 = p2(M-1);
            Parameter(1) = p1; Parameter(2) = p2;
        case 'MOEAD-DE'
            H = [99 13  7  5  4  0  3  0  2];
            H = H(M-1);
            delta = 0.9;
            nr = 2;
            Parameter(1) = H; Parameter(2) = delta; Parameter(3) = nr;
        case 'MOEAD-DRA'
            H = [99 13  7  5  4  0  3  0  2];
            H = H(M-1);
            delta = 0.9;
            Parameter(1) = H; Parameter(2) = delta;
        case 'MOEAD-M2M'
            H = [10  3  2  2  2  1  1  1  1];
            H = H(M-1);
            S = 10;
            Parameter(1) = H; Parameter(2) = S;
        case 'NSGA-II'
        case 'NSGA-III'
            p1 = [99 13  7  5  4  3  3  2  3];
            p2 = [ 0  0  0  0  1  2  2  2  2];
            p1 = p1(M-1);
            p2 = p2(M-1);
            Parameter(1) = p1; Parameter(2) = p2;
        case 'NSLS'
            miu   = 0.5;
            delta = 0.1;
            Parameter(1) = miu; Parameter(2) = delta;
        case 'PESA-II'
            div = 10;
            Parameter(1) = div;
        case 'PICEA-g'
            NGoal = M*100;
            Parameter(1) = NGoal;
        case 'SPEA2'
        case 'SPEA2-SDE'
        case 'null'
		case 'DHRS-MOEAD'
			iPara1 = [99 13 7 5 4 3 3 2 3];
			iPara2 = [ 0  0 0 0 1 2 2 2 2];
			iPara1 = iPara1(M-1);
			iPara2 = iPara2(M-1);
			Parameter(1) = iPara1; Parameter(2) = iPara2;
			iBeta = 2;
			iGamma = 20;
			Parameter(3) = iBeta; Parameter(4) = iGamma;
        otherwise
            error(['�㷨',Algorithm,'������.']);
    end
end

