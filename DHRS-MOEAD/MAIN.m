%KC,Tan----DHRS-MOEAD
function MAIN(strProblems,iM,iRun)
	clc;format compact;tic;
	%参数设置
	[iGenerations,iN,iPara1,iPara2,iBeta,iGamma] = P_settings('DHRS-MOEAD',strProblems,iM);
%--------------------------------------------------------------------
%--------------------------------------------------------------------
	%初始化向量
	iEvaluations = iGenerations * iN;
	[iN,fvectWeight] = F_weight(iPara1,iPara2,iM);
	fvectWeight(fvectWeight == 0) = 0.000001;
	iT = 20;
	iGenerations = floor(iEvaluations/iN);

	%邻居判断
	iB = zeros(iN);
	for i = 1 : iN
		for j = i : iN
			iB(i,j) = norm(fvectWeight(i,:) - fvectWeight(j,:));
			iB(j,i) = iB(i,j);
		end
	end
	[~,iB] = sort(iB,2);
	iB = iB(:,1:iT);

	%初始化种群
	[fvectPopulation,ivectBoundary,strCoding] = P_objective('init',strProblems,iM,iN);
	fvectFunctionValue = P_objective('value',strProblems,iM,fvectPopulation);
	barrXType = randi([0 1],iN,1);
	farrZmin = min(fvectFunctionValue);
	iarrCountR = zeros(iN,1);
	
	for iGene = 1 : iGenerations
		%
		fvectReParent = [];%Reference Parent Set
		fvectReOffspr = [];%Reference Offspring Set
		for i = 1 : iN
			% Using the Hybrid Recombination Strategy to product an offspring
			if iarrCountR(i) > iBeta
				barrXType(i) = mod(barrXType(i)+1,2);
			end
			farrOffspr = fnReproduction(fvectPopulation(iB(i,:)',:),barrXType(i),ivectBoundary,strCoding);
			farrOffValue = P_objective('value',strProblems,iM,farrOffspr);
			farrZmin = min(farrZmin,farrOffValue);

			for j = 1 : iT
				fFitTchPa = max(abs(fvectFunctionValue(iB(i,j),:) - farrZmin) .* fvectWeight(iB(i,j),:));
				fFitTchOff = max(abs(farrOffValue - farrZmin) .* fvectWeight(iB(i,j),:));
				if fFitTchOff < fFitTchPa
					fRDL = fnComputMRDL(fvectReParent,fvectReOffspr,fvectFunctionValue(iB(i,j),:),farrOffValue);
					iarrCountR(i) = iarrCountR(i) + 1;
					if fRDL < iGamma
						[~,iMinDisNo] = min(sum((repmat(farrOffValue,iT,1) - fvectFunctionValue(iB(i,:)',:)).^2,2));
						fvectReParent = [fvectReParent;fvectFunctionValue(iB(i,iMinDisNo),:)];
						fvectReOffspr = [fvectReOffspr;farrOffValue];
						fvectPopulation(iB(i,j),:) = farrOffspr;
						fvectFunctionValue (iB(i,j),:) = farrOffValue;
						iarrCountR(i) = 0;
					end
				end
			end
		end
		clc;fprintf('MOEA/D,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(iRun),strProblems,num2str(iM),num2str((iGene/iGenerations*100)),num2str((toc)));
	end
	P_output(fvectPopulation,toc,'DHRS-MOEAD',strProblems,iM,iRun);
end

function farrOffspr = fnReproduction(fvectPopulation,bXType,ivectBoundary,strCoding)
%混合杂交变异
	iT = size(fvectPopulation,1);
	switch bXType
	case 0
		iPareNo = randperm(iT);
		farrOffspr = P_generator(fvectPopulation(iPareNo([1;2]),:),ivectBoundary,strCoding,1);
	case 1
		iPareNo = randperm(iT);
		farrOffspr = F_generator(fvectPopulation(iPareNo(1),:),fvectPopulation(iPareNo(2),:),fvectPopulation(iPareNo(3),:),ivectBoundary);
	otherwise
	end
end


function Offspring = F_generator(r1,r2,r3,Boundary)
%微分进化,变异产生一个子代
    
    D = length(r1);
    MaxValue = Boundary(1,:);
    MinValue = Boundary(2,:);

    %微分进化参数
    CR = 1;         %控制参数
    F =0.5;       	%控制参数
    ProM = 1/D;     %变异概率
    DisM = 20;     	%变异参数
    
    %微分进化
    Offspring = r1;
    Temp = rand(1,D)<=CR;
    Offspring(Temp) = Offspring(Temp)+F.*(r2(Temp)-r3(Temp));
    
    %多项式变异
    k = rand(1,D);
    miu = rand(1,D);
    Temp = (k<=ProM & miu<0.5);
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
    Temp = (k<=ProM & miu>=0.5);
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));        
    
    %越界处理
    Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
    Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
end

function fRDL = fnComputMRDL(fvectReParent,fvectReOffspr,farrPa,farrOff)
%计算最大RDL
	if isempty(fvectReParent)
		fRDL = 0;
		return ;
	end
	iS = size(fvectReParent,1);
	fMaxRDL = 0;
	for i = 1 : iS
		fTmpRDL = fnCalculateTri(fvectReParent(i,:),fvectReOffspr(i,:),farrPa)/fnCalculateTri(fvectReParent(i,:),fvectReOffspr(i,:),farrOff);
		if fTmpRDL > fMaxRDL
			fMaxRDL = fTmpRDL;
		end
	end
	fRDL = fMaxRDL;
end

function fTriArea = fnCalculateTri(farrA,farrB,farrC)
%计算三角形面积
	iA = norm(farrA - farrB);
	iB = norm(farrB - farrC);
	iC = norm(farrC - farrA);
	iP = (iA + iB + iC) / 2;
	fTriArea = sqrt(iP * (iP - iA) * (iP - iB) * (iP - iC));
end
