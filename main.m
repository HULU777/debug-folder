%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Alp Sayin - alpsayin[at]alpsayin[dot]com - https://alpsayin.com
%   Matlab Genetic Algorithm
%   Spring 2012
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAN-UP
clear;warning('off');
% clc;close all;
tic
had = 0;
square = 1;
Q = 0;
% parameter= input(BIBD7(),had); %  = [7,7,3,9.5277];
parameter = [8 8 3 4];
My = GA_Opt(); % default values of the option field
My.numberOfGenerations = 10;  % 250
My.Z = parameter(1);
My.B = parameter(2);
My.M = parameter(3);
My.populationSize = 300; % 150 ;
My.crossoverProbability = 0.1 ;
My.mutationProbability = 0.09; % 0.0625;
My.fitness = 'Q';
My.scale = 'exp';
My.tournamentSize = 2;  %10
My.numberOfReplications = 0;  %2
My.C = 100;
My.c = 2;
My.coding = 'Hadamard';



%% RUN GENERATIONS
for iGeneration = 1: numberOfGenerations
    
    %% FIND MAXIMUM FITNESS OF POPULATION
    square = calculateS_mind(populationsquare,B,Z,had);
    Q = calculateQ_mind(populationQ,B,Z,had);
    f = figure('visible', 'off');
%     figure;
    subplot(2,1,1);
    plot(sort(square,'descend'),'r*'); hold on;
    plot(sort(Q,'descend'),'g*');
    
    scaledsquare = scalingfitness(square);
    scaledQ = scalingfitness(Q);
    subplot(2,1,2);
    plot(sort(scaledsquare,'descend'),'ro');hold on;
    plot(sort(scaledQ,'descend'),'go');
    titlename = ['BIBD8_',num2str(iGeneration),'.png'];
    saveas(f,titlename);
%     saveas(gcf,'BIBD7_%d.png',iGeneration);
    % [maximumFitness, bestIndividualIndex] = max(fitness);
    [maximumS, bestIndexS] = max(square);
    [maximumQ, bestIndexQ] = max(Q);
    xBestsquare = populationsquare(bestIndexS,:);
    xBestQ = populationQ(bestIndexS,:);
    
    [~,distancesquare] = calculateD(reshape(xBestsquare,B,Z),had);
    [~,distanceQ] = calculateD(reshape(xBestQ,B,Z),had);
    mind = min(min(distance));
    check = CHECK(mppm(4),distancesquare,distanceQ);  %(par1,par2,par)
    if check == 1
        fprintf('Maximum dmin: %d\n',par1,par2);
        fprintf('Best solution 1: \n');
        xBestclsquare = reshape(xBestsquare,B,Z);  % a column is a codeword
        disp(num2str(xBestclsquare));
        fprintf('Best solution 2: \n');
        xBestclsquare = reshape(xBestQ,B,Z);  % a column is a codeword
        disp(num2str(xBestclQ));
        break; 
    elseif check == 2
        fprintf('Maximum dmin: %d\n',par1);
        fprintf('Best solution 1: \n');
        xBestclsquare = reshape(xBestsquare,B,Z);  % a column is a codeword
        disp(num2str(xBestclsquare));
        break; 
    elseif check ==3
        fprintf('Maximum dmin: %d\n',par2);
        fprintf('Best solution 2: \n');
        xBestclQ = reshape(xBestQ,B,Z);  % a column is a codeword
        disp(num2str(xBestQ));
        break; 
    end
    
    populationsquare = operations(populationsquare,tournamentSize,scaledsquare,crossoverProbability,mutationProbability,bestIndexS,numberOfReplications,B,Z,M);
    populationQ = operations(populationQ,tournamentSize,scaledQ,crossoverProbability,mutationProbability,bestIndexQ,numberOfReplications,B,Z,M);

    
end

% Print out
% fprintf('Maximum Fitness: %d\n',maximumFitness);
% fprintf('Best Solution: %d\n',xBest); 
toc

function population = operations(population,tournamentSize,scaled,crossoverProbability,mutationProbability,bestIndividualIndex,numberOfReplications,B,Z,M)
        %% COPY POPULATION
        newPopulation = population;

        %% NEW GENERATION
        for i = 1:tournamentSize:size(population,1)
            %% TOURNAMENT SELECTION
    %         i1 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);  % ?????????
    %         i2 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
    %         chromosome1 = population(i1,:);  % 
    %         chromosome2 = population(i2,:);
            %% SCALED FITNESS SELECTION
            i1= ScaleSelect(scaled);
            i2= ScaleSelect(scaled);
            chromosome1 = population(i1,:); 
            chromosome2 = population(i2,:);

            %% CROSS-OVER
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


        %% MUTATE
        newPopulation = Mutate(newPopulation, mutationProbability,B,M,Z);
%         fprintf('iteration time: %d\n',iGeneration);

        %% check M in each section
    %     newPopulation = mutation_restore_rate(newPopulation,M,B,Z);


        %% PRESERVATION OF PREVIOUS BEST SOLUTION
        bestChromosome = population(bestIndividualIndex,:);
        newPopulation = InsertBestIndividual(newPopulation, bestChromosome, numberOfReplications);

        %% COPY THE NEW POPULATION ONTO CURRENT POPULATION
        population = newPopulation;
end

function check = CHECK(par1,par2,par)
        check = 0;
        if (par1>=par) && (par2>=par) 
            check = 1;
        end
        if (par1>=par) 
            disp('par1 finish');
            check = 2;
        end
        if (par2>=par)
            disp('par2 finish');
            check = 3;
        end
    end    