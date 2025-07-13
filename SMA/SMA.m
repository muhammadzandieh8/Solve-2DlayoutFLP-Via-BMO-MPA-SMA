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
        AllFitness(i)= Fitness(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);    
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
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
[val,idx] =sort(AllFitness);
ShowBestAnswer = 1;
temp = ShowBestAnswer;
for x=1:ShowBestAnswer
    temp(x,:)= X(idx(x),:);
    tempval(x) = val(x);
end
elapsed_time=toc;
fprintf('SMA Finished %f Seconds. \n',elapsed_time);
end

