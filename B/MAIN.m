function MAIN(Problem,M,Run)
clc;format compact;tic;
%-----------------------------------------------------------------------------
%参数设定
	[Generations,N] = P_settings('A',Problem,M);
%-----------------------------------------------------------------------------
%算法开始
	%初始化种群
	[Population,Boundary,Coding] = P_objective('init',Problem,M,N);
	%Version Alpha		
if 0
	for Gene = 1 : Generations
		if mod(N,2) == 0
			MatingPool = Population;    %交配池选择
		else
			MatingPool = [Population;Population(randi([1 N],1,1),:)];
		end
		Offspring = P_generator(MatingPool,Boundary,Coding,N);
		Population = [Population;Offspring];
		FunctionValue = P_objective('value',Problem,M,Population);
		CosDistance = ones(2*N);
		for i = 1 : 2*N
			for j = i+1 : 2*N
				CosDistance(i,j) = FunctionValue(i,:) * FunctionValue(j,:)'/(norm(FunctionValue(i,:))*norm(FunctionValue(j,:)));
				CosDistance(j,i) = CosDistance(i,j);
			end
		end
		MaxAngle = acosd(min(min(CosDistance)));
		Gamma = cosd(MaxAngle/(N/M));     %空间范围确定
		
		PopSortNo = zeros(2*N,1);
		for i = 1 : N
			Candidate = find(CosDistance(i,:)>Gamma);
			if isempty(Candidate)
				PopSortNo(i) = PopSortNo(i) + 1;
			else
				[~,TmpNo] = sort(sum(FunctionValue(Candidate,:).^2,2));
				[~,TNo] = sort(TmpNo);
				for  j = 1 : length(Candidate)
					if PopSortNo(Candidate(j)) == 0
						PopSortNo(Candidate(j)) = TNo(j);
					else
						PopSortNo(Candidate(j)) = (PopSortNo(Candidate(j)) + TNo(j)) / 2;
					end
				end
			end
		end

		[~,No] = sort(PopSortNo);

		Population = Population(No(1:N),:);
		FunctionValue = P_objective('value',Problem,M,Population);
		clc;fprintf('B,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
	end
end
	%Version Beta
if 0
	for Gene = 1 : Generations
		if mod(N,2) == 0
			MatingPool = Population;    %交配池选择
		else
			MatingPool = [Population;Population(randi([1 N],1,1),:)];
		end
		Offspring = P_generator([MatingPool],Boundary,Coding,N);
		Population = [Population;Offspring];
		FunctionValue = P_objective('value',Problem,M,Population);
		CosDistance = ones(N);
		for i = 1 : 2*N
			for j = i+1 : 2*N
				CosDistance(i,j) = FunctionValue(i,:) * FunctionValue(j,:)'/(norm(FunctionValue(i,:))*norm(FunctionValue(j,:)));
				CosDistance(j,i) = CosDistance(i,j);
			end
		end
		MaxAngle = acosd(min(min(CosDistance)));
		Gamma = cosd(MaxAngle/(N/M));     %空间范围确定
		
		PopSortNo = zeros(2*N,1);
		FrontNo = 1;
		while ~all(PopSortNo)
			for i = 1 : 2 * N
				if ~PopSortNo(i)
					Candidate = find(CosDistance(i,:)>Gamma);
					if isempty(Candidate)
						PopSortNo(i) = FrontNo; 
                    else
                        Candidate = [Candidate,i];
						[~,SortNo] = sort(sum(FunctionValue(Candidate,:).^2,2));
						for j = 1 : length(SortNo)
							if ~PopSortNo(Candidate(SortNo(j)))
								PopSortNo(Candidate(SortNo(j))) = FrontNo;
								break;
							end
						end
					end
				end
            end
            FrontNo = FrontNo + 1;
		end

		[~,No] = sort(PopSortNo);

		Population = Population(No(1:N),:);
		FunctionValue = P_objective('value',Problem,M,Population);
		clc;fprintf('B,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
	end
end

	%Version Gamma
if 0 
	FunctionValue = P_objective('value',Problem,M,Population);
	IndividualNeighbor = fnCaculateNeighbor(FunctionValue);
	for Gene = 1 : Generations
		for i = 1 : N
			MatingPool = Population;
			Offspring = P_generator(MatingPool,Boundary,Coding,1);
			OffValue = P_objective('value',Problem,M,Offspring);
			ReplaceNo = fnCalculateMinDistance(OffValue,FunctionValue);
			%Selection Strategy
			OffForce = zeros(1,M);
			PareForce = zeros(1,M);
			for j = 1 : M
				OffForce = OffForce + (OffValue - reshape(IndividualNeighbor(ReplaceNo,j,:),1,2))/(norm(OffValue-reshape(IndividualNeighbor(ReplaceNo,j,:),1,2)).^3);
				PareForce = PareForce + (FunctionValue(ReplaceNo,:) - reshape(IndividualNeighbor(ReplaceNo,j,:),1,2))/(norm(FunctionValue(ReplaceNo,:)-reshape(IndividualNeighbor(ReplaceNo,j,:),1,2)).^3);
			end
			OffRadForce =  (sum(OffForce.*ones(1,M),2)/(norm(OffForce)));
			PareRadForce =  (sum(PareForce.*ones(1,M),2)/(norm(PareForce)));
			if OffRadForce < PareRadForce
				Population(ReplaceNo,:) = Offspring;
				FunctionValue(ReplaceNo,:) = OffValue;
			end
		end
		cla;
		P_draw(FunctionValue);
		pause(1);
		clc;fprintf('B,第%2s轮,%5s问题,第%2s维,已完成%4s%%,耗时%5s秒\n',num2str(Run),Problem,num2str(M),num2str(round2(Gene/Generations*100,-1)),num2str(round2(toc,-2)));
	end
end
	%Version Delta
	if 1 
		
	end
    P_output(Population,toc,'B',Problem,M,Run);
	
end

%Calculate the Neighbors of each Point
function IndividualNeighbor = fnCaculateNeighbor(FunctionValue)
	[N,M] = size(FunctionValue);
	IndividualNeighbor = zeros(N,M,M);
	EyeMat = eye(M);
	for i = 1 : M
		[~,RankNo] = sort(FunctionValue(:,i));
		for j = 1 : N-1
			IndividualNeighbor(RankNo(j),i,:) = FunctionValue(RankNo(j+1),:);
		end
		IndividualNeighbor(RankNo(N),i,:) = EyeMat(i,:);
	end
end

%Calculate the Min Distance between Offspring and Parent
function ReplaceNo = fnCalculateMinDistance(OffValue,FunctionValue)
	[N,M] = size(FunctionValue);
	Distance = sqrt(sum((FunctionValue-repmat(OffValue,N,1)).^2,2));
	[~,ReplaceNo] = min(Distance);
end
