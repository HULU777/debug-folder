% M = 2;
% B=5;
% Z=5;
% psize = 4;
% population = initialized(psize,M,B,Z);

function population = initialize(psize,M,B,Z)  
% psize: population size
    A = rand(B,psize*Z);
    % find the column max values and location (using linear indexing) where they occur
    for m = 1:M
        [~,I] = max(A,[],1,'linear');
        % assign the max values at the locations where they occur
        A(I) = 0;
    end
    population = (A==0);
    population = reshape(population,B*Z,psize);
    population = population.';
end