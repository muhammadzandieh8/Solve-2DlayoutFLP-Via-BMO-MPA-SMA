clear all 
close all
clc

addpath('Base')
addpath('BMO')
addpath('SMA')

run = 1; % 25
Max_iteration = 2;
PopSize =2;
filename = 'result';
functionsNumber = 4;
ShowBestAnswer = 3;

solution = zeros(functionsNumber, run);
InitValues;
Answer = repmat(Chromosome(),(ShowBestAnswer*functionsNumber),MachineNumber);
currentval = 1;

optimizedPL = SMA('SMA',Max_iteration,chromosomes,PopSize, MachineNumber,LengthWorkshop,WidthWorkshop,ub,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper,f,C,ShowBestAnswer,LoC,WoC,XoC,YoC);
RESULT = BMO('BMO',Max_iteration,chromosomes,PopSize, MachineNumber,LengthWorkshop,WidthWorkshop,ub,M,L,W,Xio,Yio,Xoo,Yoo,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper,f,C,ShowBestAnswer,LoC,WoC,XoC,YoC,optimizedPL);
%DrawMap(Answer,tempval,W,L,Lo,Wo,Xo,Yo,ylower,yupper,xlower,xupper);