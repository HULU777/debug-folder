function [mindproperty,distance_eachpoint] = calculateHD(matrix)   % ,table  dmatrix  ,P
% calculate Hamming distance of a given codebook
%     P = 10^(P/10);
    [~,N] = size(matrix);
%     Cpower = sum(sum(matrix.^2));   
%     nM = matrix./ sqrt(Cpower/(P*N));    % normMatrix  n*N
    nM=matrix;
    distance_eachpoint = zeros(N,N-1);
    for n = 1:N
        points = [1:(n-1),(n+1):N];
        for j = 1: N-1
            vector1 = nM(:,n);vector2 = nM(:,points(j));
            d1 = xor(vector1 ,vector2);
            distance_eachpoint(n,j) = sum(d1);
        end
    end
    distance = reshape(distance_eachpoint,1,[]);
    distancetable = tabulate(distance);
    dmincountidx = find(distancetable(:,2));
    mindproperty =  distancetable(dmincountidx,[1,2]);
end
        
