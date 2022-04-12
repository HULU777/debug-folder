function population = operations(population,tournamentSize,scaled,selectionProbabilityMatrix,crossoverProbability,mutationProbability,bestIndividualIndex,numberOfReplications,B,Z,M)
        %% COPY POPULATION
        newPopulation = population;

        %% NEW GENERATION
        parentSize = size(population,1);
        parents = SimultaneousSelect(scaled,selectionProbabilityMatrix,parentSize);
        for i = 1:tournamentSize:size(population,1)
            %% TOURNAMENT SELECTION
    %         i1 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);  % ?????????
    %         i2 = TournamentSelect(fitness,tournamentSelectionParameter,tournamentSize);
    %         chromosome1 = population(i1,:);  % 
    %         chromosome2 = population(i2,:);
            %% SCALED FITNESS SELECTION
%             i1= ScaleSelect(scaled,selectionProbabilityMatrix);
%             i2= ScaleSelect(scaled,selectionProbabilityMatrix);
            i1 = parents(i);
            i2 = parents(i+1);
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
