
%% RUN GENERATIONS
% My = GA_Opt(); % default values of the option field
close all; clc;
My.numberOfGenerations = 10;  % 250
My.numberOfReplications = 2;
My.coding = 'None'; %'Hadamard';
My.fitnessmethod = 'Q'; %'square';
My.selectmethod = 'rankExp'; % 'exp'; % 'linear';
My.name = [My.coding,'+',My.fitnessmethod,'+',My.selectmethod];
parameter = [7,7,3,4];
My.Z = parameter(1);
My.B = parameter(2);
My.M = parameter(3);
My.optimald = parameter(4);
%     My.optimald = parameter(5);
My.populationSize = 1000; % 150 ;
My.crossoverProbability = 0.1 ;
My.mutationProbability = 0.09; % 0.0625;
My.tournamentSize = 2;  %10
My.numberOfReplications = 0;  %2
My.C = 100;
My.c = 2.7;
My.stopmethod = 'DMIN';
My.expscaleoffset = 0;
My.Qdominant = 3;
My.tournamentSelectionParameter = 0.8;
My.L = 1;
My.rankvector = [1 -2];

My1 = My;
My1.selectmethod = 'exp';
My1.fitnessmethod ='Q'; % 'Q';
My1.name = [My1.coding,'+',My1.fitnessmethod,'+',My1.selectmethod];
My1.c = 2;
My1.Qdominant =2.7273; % 7: ; 11:  3.8911;  26: 6.6667

population = initialize(My.populationSize,My.M,My.B,My.Z);
population1 = population;
population2 = population;

for iGeneration = 1: numberOfGenerations
    disp("This is generation:");
    disp(iGeneration);
    [fitness1,parents1,population1,check1] = ga1run(My1,population1);
    [fitness2,parents2,population2,check2] = ga1run(My,population2);
    edges = [0:25:1000];
%     figure; histogram(parents1,edges);hold on; histogram(parents2,edges);hold on;
    plot(sort(parents1),'r');hold on;plot(sort(parents2),'b');hold on;
%     plot2(fitness1,fitness2,parents1,parents1,iGeneration);
%     plot_3(fitness1,fitness2,scaled1,parents2,sortPselected1,sortPselected2,iGeneration);
    if check1 || check2
        break;
    end
end


function [fitness,parents,population,check] = ga1run(My,population)
    disp('Following is');
    disp(My.name);
    Z = My.Z;
    B = My.B;
    M = My.M;
    L = My.L;
    crossoverProbability = My.crossoverProbability;
    mutationProbability = My.mutationProbability; % 0.0625;
    tournamentSize = My.tournamentSize;  %10
    tournamentSelectionParameter = My.tournamentSelectionParameter;
    numberOfReplications = My.numberOfReplications;  %2
    C = My.C;
    c = My.c;
    Qdominant = My.Qdominant;
    expscaleoffset = My.expscaleoffset;
    rankvector = My.rankvector;
    if isfield(My,'optimalf')
    optimalf = My.optimalf;
    end
    if isfield(My,'optimald')
    optimald = My.optimald;
    end
    
    switch(My.coding)
        case('Hadamard'); codedpopulation = Hadamardcoding(population);
        case('None'); codedpopulation = population;
        otherwise; disp('did not find coding method'); 
    end
    
    switch (My.fitnessmethod)
        case ('Q'); [Ddistribution,fitness] = calculateQ_mind(codedpopulation,B,Z^L,C,Qdominant);
        case ('square'); [Ddistribution,fitness] = calculateS_mind(codedpopulation,B,Z^L,C,optimald);
        case ('Admin'); Ddistribution = calculate_Admin(codedpopulation,B,Z^L);  fitness = Ddistribution;
        otherwise; disp('did not find selectmethod'); 
    end
    
    switch(My.stopmethod)
        case('DMIN'); [~,Admin,bestIndex,check] = checkd(Ddistribution,optimald);
        case('fitness');[~,bestIndex,check] = checkf(fitness,optimalf);
        otherwise; disp('did not find stopmethod'); 
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
end


