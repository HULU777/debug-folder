% fitness = [2 6 8 4 5 ];
% iSelected = ScaleSelects(fitness);

function iSelected = ScaleSelect(fitness,selectionProbabilityMatrix)
%     %select 'tournamentSize' candidates for tournament
%     candidates = 1 + fix(rand(1,tournamentSize)*populationSize);
%     candidateFitnesses = fitnessValues(candidates);
%     [~, sortedIndexes] = sort(candidateFitnesses,1,'descend');
%     selectionProbabilityMatrix = fitness./sum(fitness);
%     for i = 1:length(fitness)-1
%         selectionProbabilityMatrix(i+1) = selectionProbabilityMatrix(i) + selectionProbabilityMatrix(i+1);
%     end
    r = rand;
    iSelected = find(r<selectionProbabilityMatrix);
    if isempty(iSelected)
        iSelected = length(fitness);
    else
        iSelected = iSelected(1);
    end
end