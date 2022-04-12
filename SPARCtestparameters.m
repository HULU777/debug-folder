clear; close all; clc;
% rankmatrix = [ 1 -2 0 0  0  0  0  0  0  0   0    0    0   0    0    0;
%                        1 -2 3 -4 0  0  0  0  0  0   0    0    0   0    0    0;
%                        1 -2 3 -4 5 -6  0  0  0  0   0    0    0   0    0    0;
%                        1 -2 3 -4 5 -6  7 -8  0  0   0    0    0   0    0    0;
%                        1 -2 3 -4 5 -6  7 -8  9 -10 0    0    0   0    0    0;
%                        1 -2 3 -4 5 -6  7 -8  9 -10 11 -12  0   0    0    0;
%                        1 -2 3 -4 5 -6  7 -8  9 -10 11 -12 13 -14  0    0;
%                        1 -2 3 -4 5 -6  7 -8  9 -10 11 -12 13 -14 15 -16];
% rankmatrix = [0.1;0.2;0.3;0.4];% [0.01;0.02;0.03;0.05;0.06];
% runtime = 10;
% t = ones(size(rankmatrix,1),runtime)*(-1);
% iG = ones(size(rankmatrix,1),runtime)*(-1);
% iD = ones(size(rankmatrix,1),runtime)*(-1);
% parfor i = 1:size(rankmatrix,1)
% %     pc = rankmatrix(i,:);
% %     pc = pc(pc~=0);
%     pc = rankmatrix(i);
%     for j = 1:runtime
%          [t(i,j),iG(i,j),iD(i,j)] = SPARCtestparameter(pc);
%     end
% end
 
 function [t,iGeneration,mind]=SPARCtestparameter(pc)
 %% parameters setting
% clear all; close all; clc;
% My = GA_Opt(); % default values of the option field
My.coding ='Hadamard';  % 'None'
My.fitnessmethod = 'Admin'; %'square';% 'Q';
My.selectmethod = 'rank';
% parameter = input_(BIBD7(),My.coding,My.fitnessmethod);
parameter = [7,7,3,4] ;%[Z,B,W,D] ; % [26,23,9,10];  [12,12,5,6]
My.Z = parameter(1);
My.B = parameter(2);
My.M = parameter(3);
My.optimald = parameter(4);
% My.optimalf = parameter(5);
My.populationSize = 300; % 150 ;
My.crossoverProbability = 0.3;
My.mutationProbability = 0.05; % 0.0625;
My.numberOfGenerations = 500;  % 250
My.tournamentSize = 2;  %10
My.numberOfReplications = 0;  %2
My.C = 100;
My.c = 2;
My.stopmethod = 'DMIN'; %'fitness'
My.expscaleoffset = 0;
My.Qdominant = 3;
My.numberOfReplications = 1;
My.rankvector = [1 -2];  % [1 -2 3 -4]
My.cmatrix = 0;  % HAD5_14()
My.alpha = 0.33;
My.L = 2;

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

population = initialize(populationSize,M,B,Z);

for iGeneration = 1: numberOfGenerations
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
        case ('Q'); [Ddistribution,fitness] = calculateQ_mind(codedpopulation,B,Z,C);
        case ('square'); [Ddistribution,fitness] = calculateS_mind(codedpopulation,B,Z,C,optimald);
        case ('Admin'); Ddistribution = calculate_Admin(population,B,Z);  fitness = Ddistribution;
        otherwise; disp('did not find selectmethod'); 
    end
    
    switch(My.stopmethod)
        case('DMIN'); [mind,Admin,bestIndex,check] = checkd(Ddistribution,optimald);
        case('fitness');[~,bestIndex,check] = checkf(fitness,optimalf);
    end
    
    if check
        xBest = population(bestIndex,:);
        fprintf('This is time:');
        disp(iGeneration);
        fprintf('searched optimal result: \n');
        xBestcl= reshape(xBest,B,Z);  % a column is a codeword
        disp(num2str(xBestcl));
        break;
    end
    
      switch (My.selectmethod)
        case ('exp'); [scaled,sortPselected,parents] = scalingexpSelect(fitness,c);
        case ('linear'); [scaled,sortPselected,parents] = scalinglinearSelect(fitness,c);
        case ('rank'); parents = rankSelect(fitness,c,rankvector);
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
end
    t = toc;
end

