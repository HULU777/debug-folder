%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Alp Sayin - alpsayin[at]alpsayin[dot]com - https://alpsayin.com
%   Matlab Genetic Algorithm
%   Spring 2012
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tempPopulation = Mutate(tempPopulation, mutationProbability,B,M,Z)
%% random mutation 
% indexes = rand(size(tempPopulation))<mutationProbability;                 % Index for Mutations
% tempPopulation(indexes) = tempPopulation(indexes)*-1+1;                     % Bit Flip Occurs

%% codeword mutation (0-1 swap)
tempPopulation = tempPopulation.';
pop = size(tempPopulation,2);
infoPos = find(tempPopulation==1);
infoPos = reshape(infoPos,M,pop*Z);
frozenPos = find(tempPopulation==0);
frozenPos = reshape(frozenPos,B-M,pop*Z);
% flag1 = sum(tempPopulation) ~= M;
flipinfo = randi(size(infoPos,1),1,size(infoPos,2))+(0:size(infoPos,2)-1)*size(infoPos,1);
flipfrozen = randi(size(frozenPos,1),1,size(frozenPos,2))+(0:size(frozenPos,2)-1)*size(frozenPos,1);
indexes = rand(pop*Z,1)<mutationProbability;   % corresponding codeword flip
info = infoPos(flipinfo);
frozen = frozenPos(flipfrozen);
tempPopulation(info(indexes))=0;
tempPopulation(frozen(indexes))=1;
tempPopulation = tempPopulation.';
        
% Deprecated - to be deleted in the next iteration
% nGenes= size(chromosome,2);
% mutatedChromosome = chromosome;
% for j = 1: nGenes
%     r= rand;
%     if (r < mutationProbability)
%         mutatedChromosome(j) = 1-chromosome(j);
%     end
%     
% end

