function temp = MPA(algorithmName,Max_iteration,chromosomes,PopSize, MachineNumber,LengthWorkshop,WidthWorkshop,ub,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper,f,C,ShowBestAnswer,LoC,WoC,XoC,YoC)
tic
disp('MPA is now tackling your problem...')
addpath(genpath('..'))
dim = 1;   %dimension size

lb = [0];% lower boundary 
ub = [7];    % upper boundary

% initialize position
%bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(PopSize,1);%record the fitness of all slime mold
weight = ones(PopSize,1);%fitness weight of each slime mold
Convergence_curve=zeros(1,Max_iteration);
z=0.3; % parameter
searchAgent = 1;
Xmin=repmat(ones(1,dim).*lb,searchAgent,1);
Xmax=repmat(ones(1,dim).*ub,searchAgent,1);
%disp('Initialize the set of random solutions...')
%Initialize the set of random solutions
Prey = initialization(PopSize,dim,ub,lb);;

%%%%%%%%%%%%%%%%%%%%%%%%%555
Top_predator_pos=repmat(0,1,1);
Top_predator_fit=inf; 

Convergence_curve=zeros(1,Max_iteration);
stepsize=repmat(1,PopSize,1);
fitness=inf(PopSize,1);

Iter=1;
FADs=0.2;
P=0.5;
while Iter<=Max_iteration  
disp(['MPA iterations is ' num2str(Iter)]);
     %------------------- Detecting top predator -----------------    
     for i=1:size(chromosomes,1)  
        Result = IsOverLapHappend(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);
        if Result == true
            chromosomes(i,:) = CreateCar(MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);  
        end
        %Calculate Fitness
        fitness(i) = Fitness(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);                     
     if fitness(i)<Top_predator_fit 
       Top_predator_fit=fitness(i); 
       Top_predator_pos=Prey(i,:);
     end  
    end
     %------------------- Marine Memory saving -------------------  
 if Iter==1
   fit_old=fitness;    Prey_old=Prey;
 end
     
  Inx=(fit_old<fitness);
  Indx=repmat(Inx,1,dim);
  Prey=Indx.*Prey_old+~Indx.*Prey;
  fitness=Inx.*fit_old+~Inx.*fitness;

  fit_old=fitness;    Prey_old=Prey;

     %------------------------------------------------------------   
     
 Elite=repmat(Top_predator_pos,PopSize,1);  %(Eq. 10) 
 CF=(1-Iter/Max_iteration)^(2*Iter/Max_iteration);
                             
 RL=0.05*levy(PopSize,searchAgent,1.5);   %Levy random number vector
 RB=randn(PopSize,searchAgent);          %Brownian random number vector
           
  for i=1:size(Prey,1)
     for j=1:size(Prey,2)        
       R=rand();
          %------------------ Phase 1 (Eq.12) ------------------- 
       if Iter<Max_iteration/3 
          stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));                    
          Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);          
          %--------------- Phase 2 (Eqs. 13 & 14)----------------
       elseif Iter>Max_iteration/3 && Iter<2*Max_iteration/3 
          
         if i>size(Prey,1)/2
          stepsize(i,j)=RB(i,j)*(RB(i,j)*Elite(i,j)-Prey(i,j));
          Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);    
         else
          stepsize(i,j)=RL(i,j)*(Elite(i,j)-RL(i,j)*Prey(i,j));                     
          Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);   
         end
         %----------------- Phase 3 (Eq. 15)-------------------
       else 
          stepsize(i,j)=RL(i,j)*(RL(i,j)*Elite(i,j)-Prey(i,j)); 
          Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);  
       end  
      end                                         
  end    
        
     %------------------ Detecting top predator ------------------        
 for i=1:size(Prey,1)  
        %Calculate Fitness
        fitness(i) = Fitness(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);                     
     if fitness(i)<Top_predator_fit 
       Top_predator_fit=fitness(i); 
       Top_predator_pos=Prey(i,:);
     end  
 end       
     %---------------------- Marine Memory saving ----------------   
 if Iter==0
   fit_old=fitness;    Prey_old=Prey;
 end
     
  Inx=(fit_old<fitness);
  Prey=Indx.*Prey_old+~Indx.*Prey;
  fitness=Inx.*fit_old+~Inx.*fitness;

  fit_old=fitness;    Prey_old=Prey;
     %---------- Eddy formation and FADs’ effect (Eq 16) ----------- 
                             
  if rand()<FADs
     U=rand(searchAgent,dim)<FADs;                                                                                              
     Prey=Prey+CF*((Xmin+rand(searchAgent,dim).*(Xmax-Xmin)).*U);
  else
     r=rand();  Rs=size(Prey,1);
     stepsize=(FADs*(1-r)+r)*(Prey(randperm(Rs),:)-Prey(randperm(Rs),:)); 
  end                                                      
  Iter=Iter+1;  
  Convergence_curve(Iter)=Top_predator_fit;
  % for i=1:size(Prey,1)
  %        Result = IsOverLapHappend(Prey(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);
  %       if Result == true
  %           Prey(i,:) = CreateCar(MachineNumber,LengthWorkshop,WidthWorkshop,L,W,LoC,WoC,XoC,YoC);  
  %       end
  % end
end
 for i=1:size(Prey,1)  
        fitness(i) = Fitness(chromosomes(i,:),MachineNumber,LengthWorkshop,WidthWorkshop,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,LoC,WoC,XoC,YoC,f,C);                     
 end
[val,idx] =sort(fitness);
ShowBestAnswer  = 1;
temp = repmat(searchAgent,ShowBestAnswer);
for x=1:ShowBestAnswer
    temp(x)= Prey(idx(x),:);
    tempval(x) = val(x);
end
elapsed_time=toc;
fprintf('MPA Finished %f Seconds. \n',elapsed_time);
end
