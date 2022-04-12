% c = 1.8;
% fitness = rand(10,2);
% ise = rankSelects(fitness,c);

function iSelected = rankExp(fitness,c,rankvector)
    ranklength = size(fitness,2);
    if ranklength<length(rankvector)
        rankvector = rankvector(1:ranklength);
    end
    [~,rankidx] = sortrows(fitness,rankvector);  % from bad to good  [1 -2 3 -4 5 -6]
    N = size(fitness,1);
    
    beta = c^(2/(N-1));
    alpha = (1-beta)/(beta*(1-beta^N));
    r = rand(1,N);
%     r = 0:0.01:1;
    k = log(1-r*(1-beta)/(alpha*beta))  /   log(beta);
    k = ceil(k);
%     KRANK = sort(k);
%     plot(KRANK);
    iSelected  = rankidx(k);
end