% c = 1.8;
% fitness = rand(10,2);
% ise = rankSelects(fitness,c);

function iSelected = rankLinear(fitness,c,rankvector)
    ranklength = size(fitness,2);
    if ranklength<length(rankvector)
        rankvector = rankvector(1:ranklength);
    end
    [~,rankidx] = sortrows(fitness,rankvector);  % from bad to good  [1 -2 3 -4 5 -6]
    N = size(fitness,1);
    alpha = (2*N-c*(N+1))/((N-1)*N);
    beta = 2*(c-1)/((N-1)*N);
    r = rand(1,N);
%     r = 0:0.01:1;
    k = (-(2*alpha+beta)+sqrt((2*alpha+beta)^2+8*beta*r))/(2*beta);
    k = ceil(k);
    
%     KRANK = sort(k);
%     plot(KRANK);
    iSelected  = rankidx(k);
end