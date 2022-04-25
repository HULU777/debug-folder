clear; close all; clc;

 %% parameters setting
% My = GA_Opt(); % default values of the option field
My.coding ='Hadamard';  % 'None'
My.fitnessmethod = 'Admin'; %'square';% 'Q';
My.selectmethod = 'rankExp';
% parameter = input_(BIBD7(),My.coding,My.fitnessmethod);
parameter = [16,16,2,500];%[7,7,3,40] ;%[Z,B,W,D] ; % [26,23,9,10];  [12,12,5,6]
My.Z = parameter(1);
My.B = parameter(2);
My.M = parameter(3);
My.optimald = parameter(4);
% My.optimalf = parameter(5);
My.populationSize = 300; % 150 ;
My.crossoverProbability = 0.1;
My.mutationProbability = 0.05; % 0.0625;
My.numberOfGenerations = 1000;  % 250
My.tournamentSize = 2;  %10
My.numberOfReplications = 0;  %2
My.C = 100;
My.c = 3.5;
My.stopmethod = 'DMIN'; %'fitness'
My.expscaleoffset = 0;
My.Qdominant = 3;
My.numberOfReplications = 1;
My.rankvector = [1 -2];  % [1 -2 3 -4]
My.cmatrix =HAD31_32();  % HAD4_7()  Remember change L ! eye(My.Z);  HAD5_14()
My.L = 2;
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
    
    if check || tconverge>=50 || totaltime >= 600 
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
    
    [PPM_dmin,PPM_Admin] = PPMdmin(Z,L,cmatrix,EbNo);
    disp('For PPM (same coding matrix):')
    disp('maxed_dmin:');disp(PPM_dmin);  
    disp('Admin:');disp(PPM_Admin);  
    
    semilogy(EbNo,WEF(PPM_dmin,PPM_Admin),'ro-'); hold on;
    semilogy(EbNo,WEF(mind,Admin),'b*--'); 
    title('error performance');
    legend('PPM', 'MPPM');

