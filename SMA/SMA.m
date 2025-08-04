function temp = SMA(algorithmName,Max_iteration,chromosomes,PopSize, MachineNumber,LengthWorkshop,WidthWorkshop,ub,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper,f,C,ShowBestAnswer,LoC,WoC,XoC,YoC)
tic
disp('SMA is now tackling your problem')
addpath(genpath('..'))
lb = [1];  
ub = [6];  
dim = size(ub,1);
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
%dimension size
% initialize position
%bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(PopSize,1);%record the fitness of all slime mold
weight = ones(PopSize,MachineNumber);%fitness weight of each slime mold
Convergence_curve=zeros(1,Max_iteration);
z=0.3; % parameter
N = PopSize;
X=initialization(N,dim,ub,lb);
%disp('Initialize the set of random solutions...')
%Initialize the set of random solutions

it=1;  %Number of iterations
% Main loop
while  it <= Max_iteration
    disp(['SMA iterations is ' num2str(it)]);
    %sort the fitness
    for i=1:PopSize
        Result = IsOverLapHappend(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);
        if Result == true
            chromosomes(i,:) = CreateCar(MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);  
        end
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

         % Selection (barnacle find the mating by using its penis to neighbor
         % using similar with DE
         k1 = randperm(PopSize);
         k2 = randperm(PopSize);
         
         k1x= chromosomes(k1,:);
         k2x= chromosomes(k2,:);
         pl = X(i,:);
	     % find the barnacles less than Adjusted PL
         select =[k1' k2']; 
         LessThanPLs = abs(select(:,1)-select(:,2));
         OverThanPLs = find((LessThanPLs)>pl) ;% if more than Adjusted PL the barnacle will not mating
         SamePLs = find((LessThanPLs)==0); % if 0 means self mating or spermcast
         
         Barnaclesoffspring=repmat(Chromosome(),PopSize,MachineNumber); % S_i = pop x var
         Dad_Barnacles=zeros(size(chromosomes,1),dim);
         Mom_Barnacles=zeros(size(chromosomes,1),dim);
         
         for kk=1:size(LessThanPLs,1) % 1 : PopSize  
             for kkk=1:MachineNumber
    
			     p=randn(); %p eqn 11
			     Dad_Barnacles = repmat(Chromosome(),1,MachineNumber);
			     Mom_Barnacles = repmat(Chromosome(),1,MachineNumber);
             
                 Dad_Barnacles(kk,kkk).X=p*k1x(kk,kkk).X;
                 Dad_Barnacles(kk,kkk).Y=p*k1x(kk,kkk).Y;
                 Dad_Barnacles(kk,kkk).Orientation=p*k1x(kk,kkk).Orientation;
                              
                 Mom_Barnacles(kk,kkk).X=(1-p)*k2x(kk,kkk).X;           %q=1-p  eqn 12
                 Mom_Barnacles(kk,kkk).Y=(1-p)*k2x(kk,kkk).Y;           %q=1-p  eqn 12
                 Mom_Barnacles(kk,kkk).Orientation=(1-p)*k2x(kk,kkk).Orientation; %q=1-p  eqn 12
			     % Generate new offspring 
                 Barnaclesoffspring(kk,kkk).X=XYCal((Dad_Barnacles(kk,kkk).X+Mom_Barnacles(kk,kkk).X),LengthWorkshop);
                 Barnaclesoffspring(kk,kkk).Y=XYCal((Dad_Barnacles(kk,kkk).Y+Mom_Barnacles(kk,kkk).Y),WidthWorkshop);
                 Barnaclesoffspring(kk,kkk).Orientation=OrientationCal((Dad_Barnacles(kk,kkk).Orientation+Mom_Barnacles(kk,kkk).Orientation));
             end
         end
         chromosomes(i,:) = Barnaclesoffspring(i,:);

         if OverThanPLs~=0
             for k=1:size(OverThanPLs,1)
                 temp3 = repmat(Chromosome(),1,MachineNumber);
                 temp3(k,:)=k2x(OverThanPLs(k),:);
                 for Z=1:MachineNumber             
                     temp3(k,Z).X=XYCal((rand()*temp3(k,Z).X),LengthWorkshop);
                     temp3(k,Z).Y=XYCal((rand()*temp3(k,Z).Y),WidthWorkshop);
                     temp3(k,Z).Orientation=OrientationCal(rand()*temp3(k,Z).Orientation);
				     % Generate new offspring
                     Barnaclesoffspring(OverThanPLs(k),:)=temp3(k,:);
                     chromosomes(k,:) = Barnaclesoffspring(OverThanPLs(k),:);
                 end
             end
         end
        lastfitness = AllFitness(i);
        ddd = Barnaclesoffspring(i,:);
        AllFitness(i)= Fitness(Barnaclesoffspring(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);    
        % fprintf('SMA Fitness %f  \n',AllFitness(i));

    end    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6
    worstFitness = SmellOrder(PopSize);
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
    %calculate the fitness weight of each slime mold
    for i=1:PopSize
        for j=1:size(PopSize,2)
            if i<=(PopSize/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end    
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(it/Max_iteration)+1);   %Eq.(2.4)
    b = 1-it/Max_iteration;
    % Update the Position of search agents
    for i=1:PopSize
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
            vb = unifrnd(-a,a,1,MachineNumber);  %Eq.(2.3)
            vc = unifrnd(-b,b,1,MachineNumber);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                if r<p    %Eq.(2.1)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
    it=it+1;

end
%%
% [val,idx] =sort(AllFitness);
% ShowBestAnswer = 1;
% temp = ShowBestAnswer;
% for x=1:ShowBestAnswer
%     temp(x,:)= X(val(x),:);
%     tempval(x) = val(x);
% end

% حذف مقادیر صفر از AllFitness و X مرتبط با اون‌ها
nonZeroIdx = AllFitness ~= 0;
AllFitnessClean = AllFitness(nonZeroIdx);
XClean = X(nonZeroIdx, :);

% مرتب‌سازی مقادیر غیرصفر
[val, idx] = sort(AllFitnessClean);

% نمایش بهترین پاسخ (کوچکترین مقدار غیرصفر)
ShowBestAnswer = 1;

temp = [];
tempval = [];

for x = 1:ShowBestAnswer
    if x <= length(val)  % برای جلوگیری از خطا اگر هیچ مقدار غیرصفری نباشه
        temp(x, :) = XClean(idx(x), :);
        tempval(x) = val(x);
    end
end



elapsed_time=toc;
fprintf('SMA Finished %f Seconds. \n',elapsed_time);
end

