% fitness = [2 6 8 4 5 ];
% iSelected = ScaleSelects(fitness);

function [iSelected,sortPselected] = SimultaneousscaleSelect(fitness,parentSize)
%     %select 'tournamentSize' candidates for tournament
%     candidates = 1 + fix(rand(1,tournamentSize)*populationSize);
%     candidateFitnesses = fitnessValues(candidates);
%     [~, sortedIndexes] = sort(candidateFitnesses,1,'descend');
%     selectionProbabilityMatrix = fitness./sum(fitness);
%     for i = 1:length(fitness)-1
%         selectionProbabilityMatrix(i+1) = selectionProbabilityMatrix(i) + selectionProbabilityMatrix(i+1);
%     end
    Pselected = fitness./sum(fitness);
    sortPselected = sort(Pselected);
    for i = 1:length(fitness)-1
        Pselected(i+1) = Pselected(i) + Pselected(i+1);
    end
    
    r = rand(1,parentSize);
    iSelected = ones(1,parentSize);
    for i = 1:parentSize
        smallers = find(r(i)<Pselected,1);
        if isempty(smallers)
            iSelected(i) = length(fitness);
        else
            iSelected(i) = smallers(1);
        end
    end
    %     ranki = sort(iSelected);
    %     plot(ranki)
end