%% test effect of parameters
clear;close all;clc;
run=5;
% pc = [0.1 0.5];
% pc = [1 2 3 4 5 6 7 8];  % number of replications
% pc = [500 1000 1500 2000]; % populationsize
pc = [0.01  0.09 0.2 0.5]; % iteration time
a = length(pc);
t1 = zeros(a,run);
iG1 = zeros(a,run);
f1  = zeros(a,run);
parfor j = 1: length(pc)
    for i = 1:run
        [t1(j,i),iG1(j,i),f1(j,i)] = timesga(pc(j));
    end
end



function [t,iGeneration,maximumFitness] = timesga(pc)
%% CLEAN-UP

warning('off');
tic

%% PARAMETERS
had = 0;
% mppm= input(BIBD7(),had); %  [7,7,3,9.5277];  % [11,11,5,4.2858]
mppm = [8, 8, 3, 4];
populationSize = 300; % 150 ;
crossoverProbability = 0.1; %0.09 ;
mutationProbability = 0.01; % 0.0625;
tournamentSelectionParameter = 0.5;
numberOfGenerations = 50;  % 250
tournamentSize = 2;  %10
numberOfReplications =8;  %2

Z = mppm(1);
n = 7;
R = log2(Z)/n;
B = mppm(2);
M = mppm(3);


verbose = 0;% true;
draw_plots = 0; % true;
% numberOfGenes = 40;
% variableRange = 10.0;
% numberOfVariables = 2;
% UNLESS THE FITNESS FUNCTION IS REALLY DIFFICULT TO COMPUTE, IT'S FASTER
% NOT TO USE PARALLEL COMPUTATION
% runparallel = false; 
% best773 = [1 1 0 1 0 0 0;
%     0 1 1 0 1 0 0;
%     0 0 1 1 0 1 0;
%     0 0 0 1 1 0 1;
%     1 0 0 0 1 1 0;
%     0 1 0 0 0 1 1 ;
%     1 0 1 0 0 0 1] ;
% best = best773;

%% VARIABLES
% fitness = zeros(populationSize, 1);

%% PLOTTING SETUP
if draw_plots
    fitnessFigureHandle = figure;
    hold on;
    set(fitnessFigureHandle,'Position',[50,50,500,200]);
    set(fitnessFigureHandle,'DoubleBuffer','on');
    axis([1 numberOfGenerations -variableRange variableRange]);
    bestPlotHandle = plot(1:numberOfGenerations, zeros(1,numberOfGenerations));
    textHandle = text(30,2.6, sprintf('best: %4.3f', 0.0));
    hold off;
    drawnow;

    surfaceFigureHandle= figure;
    hold on;
    set(surfaceFigureHandle,'DoubleBuffer','on');
    delta=0.1;
    limit = fix(2*variableRange/delta)+1 ;
    [xValues, yValues] = meshgrid(-variableRange: delta:variableRange,-variableRange: delta:variableRange);
    zValues= zeros(limit,limit);
    for j = 1: limit
        for k = 1: limit
            zValues(j,k) = EvaluateIndividual([xValues(j,k) yValues(j,k)]);
        end
    end
    surfl(xValues,yValues,zValues)
    colormap gray;
    shading interp;
    view ([-7 -9 10]);
    decodedPopulation = zeros(populationSize,numberOfVariables);
    populationPlotHandle = plot3(decodedPopulation(:,1),decodedPopulation(:,2),fitness(:),'kp');
    hold off;
    drawnow;
end


%% INITIATE POPULATION  (?? initialize)
% population = InitializePopulation(populationSize, numberOfGenes) ;   
population = initialize(populationSize,M,B,Z);


%% RUN GENERATIONS
for iGeneration = 1: numberOfGenerations
    
    %% FIND MAXIMUM FITNESS OF POPULATION

%     decodedPopulation = DecodePopulation(population, numberOfVariables, variableRange);
    [~,fitness] = calculate_mind(population,B,Z,had); % EvaluatePopulation(decodedPopulation, runparallel);
    scaled = scalingfitness(fitness);
    [maximumFitness, bestIndividualIndex] = max(fitness);
    xBest = population(bestIndividualIndex,:); %decodedPopulation(bestIndividualIndex,:);
    [~,distance] = calculateD(reshape(xBest,B,Z),had);
    mind = min(min(distance));
    if (mind>=mppm(4))
        fprintf('Succeed! Maximum dmin: %d\n',mind);
        fprintf('Best solution: \n');
        xBestcl = reshape(xBest,B,Z);  % a column is a codeword
        disp(num2str(xBestcl));
        break; 
    end
%     if (maximumFitness>=mppm(4))
%         xBestcl = reshape(xBest,B,Z); 
%         disp(num2str(xBestcl));
%         break; end
%   % Deprecated - to be deleted in the next iteration
%     maximumFitness = 0.0;
%     for i = 1: populationSize
%         chromosome = population(i,:);
%         x = DecodeChromosome(chromosome, numberOfVariables, variableRange) ;
%         decodedPopulation(i,:)= x;
%         fitness(i) = EvaluateIndividual(x);
%         if ( fitness(i)> maximumFitness)
%             maximumFitness = fitness(i);
%             bestIndividualIndex = i;
%             xBest=x ;
%         end
%     end

    % Print out
    if verbose
        fprintf('Maximum Fitness: %d\n',maximumFitness);
        fprintf('Best solution: \n');
        xBestcl = reshape(xBest,B,Z);  % a column is a codeword
        xBestrow = xBestcl.'; % a row is a codeword
        disp(num2str(xBestcl));
    end
    
    %% COPY POPULATION
    newPopulation = population;

    %% NEW GENERATION
    for i = 1:tournamentSize:populationSize
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
    fprintf('iteration time: %d\n',iGeneration);
    
    %% check M in each section
    newPopulation = mutation_restore_rate(newPopulation,M,B,Z);
    
    
    %% PRESERVATION OF PREVIOUS BEST SOLUTION
    bestChromosome = population(bestIndividualIndex,:);
    newPopulation = InsertBestIndividual(newPopulation, bestChromosome, numberOfReplications);
        
    %% COPY THE NEW POPULATION ONTO CURRENT POPULATION
    population = newPopulation;
    
    
    %% PLOT CURRENT SITUATION
    if draw_plots
        plotvector = get(bestPlotHandle,'YData');
        plotvector(iGeneration)= maximumFitness;
        set(bestPlotHandle,'YData',plotvector);
        set(textHandle,'String', sprintf('best: %4.3f',maximumFitness));
        set(populationPlotHandle,'XData', decodedPopulation(:,1),'YData',decodedPopulation(:,2),'ZData', fitness(:));
        drawnow;
    end
end
xBest = reshape(xBest,B,Z);
t = toc;
end

function mppm = input(optimalmppm,had)
    mppm = ones(1,4)*(-1);  % [Z,B,M,fitness]
    mppm(1:2) = size(optimalmppm);
    mppm(3) = sum(sum(optimalmppm))/mppm(1);
    disp('please be reminded each row is a codeword!!');   
    population = reshape(optimalmppm.',1,[]);
    if had ==0
        [~,mppm(4)] = calculate_mind(population,mppm(2),mppm(1),had);
    else
        mppm(4) = 20;
    end
end

function bibd = BIBD7()
bibd = [0    1    1    0    0    1    0
   0    0    0    0    1    1    1
   0    0    1    1    0    0    1
  0    1    0    1    1    0    0
   1    1    0    0    0    0    1
   1    0    0    1    0    1    0
  1    0    1    0    1    0    0];
end

function bibd = BIBD11()
 bibd = [  0    1    0    0    1    1    1    0    0     0     1
               0    0    0    0    0    0    1    1    1     1     1
               0    0    0    1    1    1    0    0    1     1     0
               0    1    1    1    0    0    0    0    1     0     1
               1    0    1    0    0    1    0    0    0     1     1
               1    0    0    1    1    0    0    1    0     0     1
               1    0    1    0    1    0    1    0    1     0     0
               0    1    1    0    1    0    0    1    0     1     0
               0    0    1    1    0    1    1    1    0     0     0
               1    1    0    1    0    0    1    0    0     1     0
               1    1    0    0    0    1    0    1    1     0     0];
end

function bibd = BIBD21()
bibd1 = [0    0    0    0    0    0    0    1    0     0     1     1     0
    0    0    0    0    0    0    1    0    0     1     0     1     0
    0    1    1    0    0    0    0    0    0     0     0     0     0
   0    0    0    1    0    0    0    0    0     0     0     0     1
    1    1    0    0    0    0    0    0    0     0     0     1     1
  0    0    1    1    1    0    0    0    0     0     0     1     0
   0    0    1    0    0    0    0    0    0     1     1     0     1
  0    0    0    1    0    0    0    0    1     0     1     0     0
    0    1    0    0    1    0    0    0    1     1     0     0     0
   0    1    0    1    0    0    1    1    0     0     0     0     0
    0    1    0    0    0    1    0    0    0     0     1     0     0
    1    0    0    1    0    1    0    0    0     1     0     0     0
    0    0    0    0    0    1    0    0    1     0     0     1     0
   0    0    0    0    0    0    0    1    0     1     0     0     0
   1    0    0    0    1    0    1    0    0     0     1     0     0
   0    0    0    0    1    0    0    0    0     0     0     0     0
   1    0    0    0    0    0    0    0    0     0     0     0     0
    1    0    1    0    0    0    0    1    1     0     0     0     0
   0    0    0    0    0    0    1    0    1     0     0     0     1
   0    0    1    0    0    1    1    0    0     0     0     0     0
   0    0    0    0    1    1    0    1    0     0     0     0     1]; 
bibd2 = [ 0     0     0     0     0     1     0     1
    1     0     0     0     1     0     0     0
    0     1     0     0     1     0     0     1
     1     1     0     0     0     1     0     0
    0     0     0     1     0     0     0     0
    0     0     1     0     0     0     0     0
    0     0     0     0     0     0     1     0
    0     0     0     1     1     0     0     0
    0     0     0     0     0     1     0     0
    0     0     0     0     0     0     1     0
     1     0     1     0     0     0     0     0
     0     0     0     0     0     0     0     1
    0     1     0     0     0     0     1     0
    0     1     1     1     0     0     0     0
   0     1     0     0     0     0     0     0
    1     0     0     1     0     0     1     1
     0     0     1     0     1     1     1     0
   1     0     0     0     0     0     0     0
    0     0     1     0     0     0     0     1
   0     0     0     1     0     1     0     0
    0     0     0     0     1     0     0     0
];
bibd = [bibd1 bibd2];
end

function bibd = BIBD31()
bibd1 = [1    0    0    1    0    1    0    0    0     1     0     0     0
    0    0    0    0    0    0    0    0    0     0     0     0     0
   1    1    1    0    0    0    0    0    0     0     0     0     0
   0    0    0    0    1    1    1    0    0     0     0     0     0
   0    0    0    1    1    0    0    0    0     0     0     0     0
   1    0    0    0    0    0    0    0    0     0     0     0     1
   0    0    1    0    0    0    1    0    0     1     0     0     0
    0    0    1    0    0    0    0    1    0     0     0     0     0
    0    0    0    0    0    1    0    0    0     0     0     0     0
    1    0    0    0    1    0    0    1    1     0     1     0     0
   0    1    0    1    0    0    0    0    0     0     0     1     1
   1    0    0    0    0    0    0    0    0     0     0     0     0
   0    1    0    0    0    0    0    0    0     0     1     0     0
    0    0    0    0    0    0    1    0    1     0     0     0     1
    0    0    0    0    1    0    0    0    0     0     0     1     0
    0    0    0    0    0    0    0    0    0     0     0     0     0
    0    0    0    0    0    0    0    0    0     0     1     0     0
    0    0    0    0    0    0    0    0    1     1     0     0     0
   0    0    0    0    0    0    0    1    0     1     0     1     0
   0    0    1    1    0    0    0    0    1     0     0     0     0
   1    0    0    0    0    0    1    0    0     0     0     1     0
   0    0    0    0    0    0    0    0    1     0     0     1     0
    0    0    0    0    0    1    0    1    0     0     0     0     1
   0    0    0    1    0    0    1    0    0     0     1     0     0
   0    1    0    0    1    0    0    0    0     1     0     0     0
   0    1    0    0    0    1    0    0    1     0     0     0     0
   0    0    1    0    1    0    0    0    0     0     0     0     1
   0    0    1    0    0    1    0    0    0     0     1     1     0
    0    1    0    0    0    0    1    1    0     0     0     0     0
    0    0    0    1    0    0    0    1    0     0     0     0     0
   0    0    0    0    0    0    0    0    0     1     1     0     1];
bibd2 = [ 0     0     0     0     0     1     0     0     0     0     0     0
     0     0     0     0     0     0     0     0     0     0     0     0
   0     0     0     0     0     0     0     0     0     0     1     1
    0     0     0     0     0     0     0     0     0     1     0     1
    0     0     0     0     0     0     0     1     1     0     1     0
     0     0     0     0     1     0     1     0     1     0     0     0
    0     0     0     0     0     0     1     1     0     0     0     0
    0     0     0     0     1     1     0     0     0     1     0     0
     0     0     1     1     1     0     0     0     0     0     1     0
    0     0     0     0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0     1     0     0
    1     0     0     1     0     0     0     1     0     1     0     0
     0     0     1     0     0     1     0     1     0     0     0     0
    0     0     0     0     0     1     0     0     0     0     1     0
     0     0     0     1     0     1     1     0     0     0     0     0
   1     1     0     0     0     1     0     0     1     0     0     1
     0     1     0     0     0     0     1     0     0     1     1     0
    0     0     1     0     0     0     0     0     1     1     0     0
    1     0     0     0     0     0     0     0     0     0     1     0
     0     1     0     1     0     0     0     0     0     0     0     0
    0     1     1     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     1     0     0     0     1
    0     1     0     0     0     0     0     1     0     0     0     0
     1     0     0     0     1     0     0     0     0     0     0     0
    0     1     0     0     1     0     0     0     0     0     0     0
    1     0     0     0     0     0     1     0     0     0     0     0
   1     0     1     0     0     0     0     0     0     0     0     0
   0     0     0     0     0     0     0     0     1     0     0     0
    0     0     0     1     0     0     0     0     1     0     0     0
    0     0     1     0     0     0     1     0     0     0     0     1
    0     0     0     1     0     0     0     0     0     0     0     1];
bibd3 = [ 0     0     0     1     0     0
     1     1     1     1     1     1
    0     0     0     0     0     1
    0     0     0     0     1     0
    0     0     1     0     0     0
    0     0     0     0     1     0
   1     0     0     0     0     0
    0     0     1     0     0     0
    1     0     0     0     0     0
     1     0     0     0     0     0
     1     0     0     0     0     0
     0     1     0     0     0     0
     0     0     0     0     1     0
     0     1     0     0     0     0
    0     0     0     0     0     1
    1     0     0     0     0     0
     0     0     0     1     0     0
    0     0     0     0     0     1
     0     0     0     0     1     0
    0     0     0     0     1     0
     0     0     1     0     0     0
    0     0     0     1     0     0
    0     0     0     0     0     1
     0     0     0     0     0     1
     0     1     0     0     0     0
     0     0     1     0     0     0
     0     0     0     1     0     0
    0     1     0     0     0     0
     0     0     0     1     0     0
     0     1     0     0     0     0
    0     0     1     0     0     0
];
bibd = [bibd1 bibd2 bibd3];

end

function cw = CW16()
cw = [1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0
1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0
1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0
1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
0 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0
0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1
0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0
0 0 1 0 0 0 0 0 0 0 1 0 0 1 0 0
0 0 0 0 1 0 0 0 0 0 0 0 1 0 1 0
0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0
0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 1 0 0 1 0 0 0 0 1 0
0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 1
0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0];
end

%% best CWC 26
function bestCWC26 = CWC26()
bestCWC26 = [1 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
             0 1 0 1 0 0 0 1 1 1 0 0 0 0 1 0 1 0 1 0 0 1 0;
             0 0 1 0 1 1 1 0 1 0 1 0 0 1 0 1 0 0 0 1 0 0 0;
             1 0 1 1 1 0 0 0 0 1 0 0 0 0 0 1 1 1 0 0 0 0 1;
             1 0 0 0 0 0 0 0 1 1 0 1 0 1 0 0 1 1 1 0 1 0 0;
             0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 0 1 0 1 1 1 1 1;
             0 0 0 1 0 1 0 0 1 1 1 0 1 0 0 1 0 0 1 0 1 0 0;
             0 1 0 0 0 0 0 1 1 1 1 1 0 1 0 1 0 0 0 0 0 0 1;
             0 0 0 0 1 1 0 0 0 0 1 0 0 0 1 1 0 1 0 0 1 1 1;
             0 0 0 0 1 0 1 0 1 0 0 1 1 0 1 1 1 0 0 0 0 0 1;
             1 0 0 0 0 1 1 0 1 1 0 0 0 0 1 0 0 0 0 1 1 1 0;
             0 0 0 0 1 0 0 0 0 1 1 0 1 1 1 0 0 1 1 1 0 0 0; 
             0 1 0 1 0 1 1 0 0 0 0 0 1 1 0 0 1 1 0 0 1 0 0;
             0 0 0 1 0 0 0 1 1 0 1 0 0 0 0 0 1 1 0 1 1 0 1;
             0 0 1 0 1 0 1 1 0 1 0 1 1 0 0 0 0 1 0 0 1 0 0;
             1 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0 1 0 1 0;
             0 0 1 0 0 0 0 0 1 0 0 0 1 1 0 1 0 1 1 0 0 1 1;
             1 1 0 0 0 1 1 0 1 0 1 0 0 0 0 0 1 0 1 0 0 0 1;
             0 0 1 1 0 1 1 0 0 1 0 1 0 1 0 0 0 0 0 0 0 1 1;
             0 0 1 0 0 1 0 1 1 0 1 1 0 0 1 0 0 1 1 0 0 0 0;
             1 0 1 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 1 0 1 1;
             0 1 1 1 0 0 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 0;
             0 1 1 0 0 1 0 1 0 1 0 0 1 0 1 0 0 0 0 1 0 0 1;
             0 0 0 1 1 0 0 1 0 0 0 1 0 1 0 1 1 0 0 0 1 1 0;
             1 0 0 1 1 0 1 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 1;
             1 1 1 0 0 0 0 0 0 0 1 0 0 1 1 1 1 0 0 0 1 0 0;
];
end

function bestCWC61 = CWC61()
bestCWC61 = [1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
                        1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0
                        0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0
                        0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 1 1 0 0
                        0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1
                        0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0
                        0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0
                        0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 1 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1
                        0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 1 0
                        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0
                        0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0
                        0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1
                        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
                        0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0
                        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0];
end        