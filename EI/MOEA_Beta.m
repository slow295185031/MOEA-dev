%测试电子受力平衡理论
function MOEA_Beta(strProblems,iM,iRun)
clc;format compact;tic;
%------------------------------------------------------------------
%参数设定
	[iGenerations,iN] = P_settings('A',strProblems,iM);
%------------------------------------------------------------------
	%初始化种群
	[farrPopulation,farrBoundary,strCoding] = P_objective('init',strProblems,iM,iN);
%	farrFunctionValue = P_objective('value',strProblems,iM,farrPopulation);
	iTestModel = 1;
	for iGene = 1 : iGenerations
		farrOffspring = P_generator(farrPopulation,farrBoundary,strCoding,iN);
		farrPopulation = [farrPopulation;farrOffspring];
		farrFunctionValue = P_objective('value',strProblems,iM,farrPopulation);
		[fAbsRad,fAbsTang] = fnCalculateVector(farrFunctionValue);
		iarrSequence = fnSortModel(fAbsRad,fAbsTang,iTestModel);
		farrPopulation = farrPopulation(find(iarrSequence<=iN/2),:);
		%This is in processing drawing
		farrFunctionValue = P_objective('value',strProblems,iM,farrPopulation);
		cla;
		P_draw(farrFunctionValue);
		pause();
	end
	P_output(farrPopulation,toc,'EI',strProblems,iM,iRun);
end


%计算目标函数值的局部受力
function [fAbsRad,fAbsTang] =  fnCalculateVector(farrFunctionValue)
%------------------------------------------------------------------
%------------------------------------------------------------------
	[iN,iM] = size(farrFunctionValue);
	farrForce = zeros(iN,iM);
	F = zeros(iN,iM);
	for iDim = 1 : iM
		[~,iValueRank] = sort(farrFunctionValue(:,iDim));
		farrForce(iValueRank(end),iDim) = 0;
		for iNo = 1 : iN-1 
			farrForce(iValueRank(iNo),iDim) = iValueRank(iNo + 1);
			F = F + farrFunctionValue()
		end
	end
end

%This is for testing the selecting Machanism
function iarrSequence = fnSortModel(fAbsRad,fAbsTang,iTestModel)
%------------------------------------------------------------------
%------------------------------------------------------------------
	switch(iTestModel)
	case 1
		fIndicator = fAbsRad + 1 * fAbsTang;
		[~,iarrSequence] = sort(fIndicator);
	case 2
		
	otherwise
	end
end
