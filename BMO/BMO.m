function temp = BMO(algorithmName,Max_iteration,chromosomes,PopSize, MachineNumber,LengthWorkshop,WidthWorkshop,ub,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper,f,C,ShowBestAnswer,LoC,WoC,XoC,YoC,optimizedPL)
tic
disp('BMO is now tackling your problem')
dimension = 3;
dim=dimension;
pl=optimizedPL;
addpath(genpath('..'))
lb = ones(1,dim);
tic

flag=0;
if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

%Initialize the population of barnacles
BarnaclesPositions = chromosomes;
BarnaclesFitness = zeros(1,PopSize);
B_Fitness  =  zeros(1,PopSize);
fitness_history=zeros(PopSize,Max_iteration);
position_history=repmat(Chromosome(),PopSize,MachineNumber);
Convergence_curve=zeros(1,Max_iteration);


%Calculate the fitness of  barnacles
for i=1:size(BarnaclesPositions,1)
    BarnaclesFitness(i)= Fitness(BarnaclesPositions(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);
    fitness_history(i)=BarnaclesFitness(i);
    position_history(i,:)=BarnaclesPositions(i,:);
end
[sorted_fitness,sorted_indexes]=sort(BarnaclesFitness);
% Find the best barnacle in the first population 
for newindex=1:PopSize
    Sorted_barnacle(newindex,:)=BarnaclesPositions(sorted_indexes(newindex),:);
end
BarnaclesPositions=Sorted_barnacle;
BarnaclesFitness=sorted_fitness;
TargetPosition=Sorted_barnacle(1,:);
TargetFitness=sorted_fitness(1);
Convergence_curve(1)=TargetFitness;
Iteration = 1;
% Main loop
 % Start from the second iteration since the first iteration was dedicated to calculating the fitness of barnacle
while Iteration<=Max_iteration
     disp(['BMO iterations is ' num2str(Iteration)]);
     Iteration = Iteration + 1;
        
     % Selection (barnacle find the mating by using its penis to neighbor
     % using similar with DE
     k1 = randperm(PopSize);
     k2 = randperm(PopSize);
     
     k1x= BarnaclesPositions(k1,:);
     k2x= BarnaclesPositions(k2,:);
     
	 % find the barnacles less than Adjusted PL
     select =[k1' k2']; 
     LessThanPLs = abs(select(:,1)-select(:,2));
     OverThanPLs = find((LessThanPLs)>pl) ;% if more than Adjusted PL the barnacle will not mating
     SamePLs = find((LessThanPLs)==0); % if 0 means self mating or spermcast
     
     Barnaclesoffspring=repmat(Chromosome(),PopSize,MachineNumber); % S_i = pop x var
     Dad_Barnacles=zeros(size(BarnaclesPositions,1),dim);
     Mom_Barnacles=zeros(size(BarnaclesPositions,1),dim);
     
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
             end
         end
     end
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BarnaclesPositions
    
    % concatenate Parent and offspring
    Barnaclescolony= [BarnaclesPositions;Barnaclesoffspring];
    
    % Relocate barnacles that go outside the search space 
    for i=1:size(Barnaclescolony,1)
        Result = IsOverLapHappend(Barnaclescolony(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);
        if Result == true
            Barnaclescolony(i,:) = CreateCar(MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);  
        end
        
         % Calculating the objective values for all barnacles
         BarnaclesFitness(i)= Fitness(Barnaclescolony(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);        %Iteration = Iteration + 1;
    end        
   
    % sorting   
    [sorted_fitness,sorted_indexes]=sort(BarnaclesFitness);
    
        for newindex=1:PopSize*2
            Sorted_barnacle(newindex,:)=Barnaclescolony(sorted_indexes(newindex),:);
        end

        %select the best top N
		%Sorted_barnacle
        BarnaclesPositions=Sorted_barnacle(1:PopSize,:);
        SortfitbestN = sorted_fitness(1:PopSize);           
        fitness_history(:,Iteration)=SortfitbestN';
        position_history(:,:,Iteration)=BarnaclesPositions;
        %Update the target        
        if SortfitbestN(1)<TargetFitness
            TargetPosition=BarnaclesPositions(1,:);
            TargetFitness=SortfitbestN(1);
        end
    Convergence_curve(Iteration)=TargetFitness;
end

if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end
bestSolution = TargetPosition;
bestFitness = TargetFitness;
iteration = Iteration;
for i=1:size(Barnaclescolony,1)
    BarnaclesFitness(i)= Fitness(Barnaclescolony(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);
end     
[val,idx] =sort(BarnaclesFitness);
temp = repmat(Chromosome(),ShowBestAnswer,MachineNumber);
for x=1:ShowBestAnswer
    temp(x,:)= Barnaclescolony(idx(x),:);
    tempval(x) = val(x);
end
elapsed_time=toc;
DrawMap(algorithmName,LengthWorkshop,WidthWorkshop,temp,tempval,W,L,LoC,WoC,XoC,YoC);
fprintf('BMO Finished %f Seconds. \n',elapsed_time);
end