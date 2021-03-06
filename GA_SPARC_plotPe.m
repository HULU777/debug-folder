clear; close all; clc;

 %% parameters setting
EbNoRange =-3:3:21 ; EbNoRange = EbNoRange.';   % 20:27 15:20
runtime = 1;
countdmin = 1;  % count ranking the first countdmin dmins
parameter = [7,7,0,2,40];  %[Z,B,W,L,D] ;
Z = parameter(1);
B = parameter(2);
L = parameter(4);
D = parameter(5);
weight = 1:floor(parameter(2)/2); weight = weight.';
% weight = 1;
mppmA = ones(size(weight,1),runtime)*(-1);
mppmD = ones(size(weight,1),runtime)*(-1);
MPPM_dmin = ones(countdmin,length(EbNoRange),length(weight))*(-1);
MPPM_Admin = ones(countdmin,length(EbNoRange),length(weight))*(-1);
Pe = ones(length(weight),length(EbNoRange))*(-1);
Pe1 = ones(1,length(EbNoRange))*(-1);
parfor i = 1:length(weight)
%     pc = rankmatrix(i,:);
%     pc = pc(pc~=0);
    %parameter(3) = weight(i);
    w = weight(i);
    for j = 1:runtime
         [xBestcl,My,Ddistribution] = GASPARC(Z,B,w,L,D);
         [MPPM_dmin(:,:,i),MPPM_Admin(:,:,i)] = Dmin_for_EbNo(My.cmatrix,xBestcl,My.L,EbNoRange,countdmin);
         Pe(i,:) = WEF(MPPM_dmin(:,:,i), MPPM_Admin(:,:,i),Z^L);
    end
end

figure;
for i= 1:length(weight)
    semilogy(EbNoRange,Pe(i,:)); hold on;
end
title('[7,7,2],L = 2; HAD13\_14');
% legend on;
% legend('PPM', 'MPPM');


 
 function [xBestcl,My,Ddistribution]=GASPARC(Z,B,w,L,D)
% My = GA_Opt(); % default values of the option field
My.coding ='Hadamard';  % 'None'
My.fitnessmethod = 'Admin'; %'square';% 'Q';
My.selectmethod = 'rankExp';
% parameter = input_(BIBD7(),My.coding,My.fitnessmethod);
% parameter = [7,7,3,40] ;%[Z,B,W,D] ; % [26,23,9,10];  [12,12,5,6]
My.Z = Z;
My.B = B;
My.M = w;
My.optimald = D;
% My.optimalf = parameter(5);
My.populationSize = 300; % 150 ;
My.crossoverProbability = 0.1;
My.mutationProbability = 0.05; % 0.0625;
My.numberOfGenerations = 200;  % 250
My.tournamentSize = 2;  %10
My.numberOfReplications = 0;  %2
My.C = 100;
My.c = 3.5;
My.stopmethod = 'DMIN'; %'fitness'
My.expscaleoffset = 0;
My.Qdominant = 3;
My.numberOfReplications = 1;
My.rankvector = [1 -2];  % [1 -2 3 -4]
My.cmatrix =HAD13_14();  % HAD4_7()  Remember change L ! eye(My.Z);  HAD5_14()
My.L = L ;
My.alpha = 1;%0.33;
My.EbNo = 1;  % (dB)


%% RUN GENERATIONS
tic
Z = My.Z;
B = My.B;
M = My.M;
crossoverProbability = My.crossoverProbability;
mutationProbability = My.mutationProbability; % 0.0625;
tournamentSize = My.tournamentSize;  %10
numberOfReplications = My.numberOfReplications;  %2
numberOfGenerations = My.numberOfGenerations;
populationSize = My.populationSize;
C = My.C;
c = My.c;
% optimalf = My.optimalf;
optimald = My.optimald;
Qdominant = My.Qdominant;
rankvector = My.rankvector;
cmatrix = My.cmatrix;
alpha = My.alpha;
L = My.L;
EbNo = My.EbNo;

population = initialize(populationSize,M,B,Z);
oldmind = 0; oldAdmin = 0; tconverge = 0; totaltime = 0;

for iGeneration = 1: numberOfGenerations
    tic
    switch(My.coding)
        case('Hadamard')
            if cmatrix == 0
               [codedpopulation,cmatrix] = Hadamardcoding(population,Z,B,alpha,L,cmatrix);
            else
               [codedpopulation,~] = Hadamardcoding(population,Z,B,alpha,L,cmatrix); 
            end
        case('None'); codedpopulation = population;
    end
    
    switch (My.fitnessmethod)
        case ('Q'); [Ddistribution,fitness] = calculateQ_mind(codedpopulation,B,Z^L,C);
        case ('square'); [Ddistribution,fitness] = calculateS_mind(codedpopulation,B,Z^L,C,optimald);
        case ('Admin'); Ddistribution = calculate_Admin(codedpopulation,B,Z^L,EbNo);  fitness = Ddistribution;
        otherwise; disp('did not find selectmethod'); 
    end
    
    switch(My.stopmethod)
        case('DMIN'); [mind,Admin,bestIndex,check] = checkd(Ddistribution,optimald);
        case('fitness');[~,bestIndex,check] = checkf(fitness,optimalf);
    end
    
    if mind == oldmind && oldAdmin == Admin
        tconverge = tconverge+1;
    else
        tconverge = 0;
    end
    oldmind = mind;
    oldAdmin = Admin;
    
    if check || tconverge>=100 || totaltime >= 600% 7200 control the time per search
        disp(check);disp(tconverge);disp(totaltime);
        xBest = population(bestIndex,:);
        fprintf('This is time:');
        disp(iGeneration);
        fprintf('searched optimal result: \n');
        xBestcl= reshape(xBest,B,Z);  % a column is a codeword
        disp(num2str(xBestcl));
%         pause;
        break;
    end
    
      switch (My.selectmethod)
        case ('exp'); [scaled,sortPselected,parents] = scalingexpSelect(fitness,c);
        case ('linear'); [scaled,sortPselected,parents] = scalinglinearSelect(fitness,c);
        case ('rankLinear'); parents = rankLinear(fitness,c,rankvector);
        case ('rankExp'); parents = rankExp(fitness,c,rankvector);
        case ('tournament'); parents = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
        otherwise; disp('did not find selectmethod'); 
%             case ('AMPseededFourier'); AMPseededFourier();
     end
    
    newPopulation = population;
    for i = 1:tournamentSize:size(population,1)
        i1 = parents(i);
        i2 = parents(i+1);
        chromosome1 = population(i1,:); 
        chromosome2 = population(i2,:);
        r = rand;
        if ( r < crossoverProbability)
            newChromosomePair = Cross(chromosome1, chromosome2,B,Z);
            newPopulation(i,:) = newChromosomePair(1,:);
            newPopulation(i+1,:) = newChromosomePair(2,:);
        else  % no cross
            newPopulation(i,:) = chromosome1;
            newPopulation(i+1,:) = chromosome2;
        end
    end
    newPopulation = Mutate(newPopulation, mutationProbability,B,M,Z);
    bestChromosome = population(bestIndex,:);
    newPopulation = InsertBestIndividual(newPopulation, bestChromosome, numberOfReplications);

    population = newPopulation;
    t1 = toc;
    totaltime = t1+totaltime;
end
    t = toc;
 end
    
