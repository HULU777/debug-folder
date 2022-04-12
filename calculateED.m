function [mindproperty,distance] = calculateED(matrix)   % ,table  dmatrix  ,P
% calculate the Euclidean mind of a given codebook
%     P = 10^(P/10);
    P = 1; % power per blocklength
    [nn,N] = size(matrix);
    mindproperty = ones(N,2)*(-1);
    
    % normalize the average power of codeword to P=n
    % thus, total codebook power is nN. 
   
    Cpower = sum(sum(matrix.^2));   
    nM = matrix./ sqrt(Cpower/(P*nn));    % normMatrix  n*N

    
%     nM=matrix;
    distance_eachpoint = zeros(N,N-1);
    for n = 1:N
        points = [1:(n-1),(n+1):N];
        for j = 1: N-1
            d1 = nM(:,n)-nM(:,points(j));
            distance_eachpoint(n,j) = sqrt(sum(d1.^2,1));
        end
        distance = reshape(distance_eachpoint,1,[]);
        dmin = min(distance);
        mindproperty(1,1) = dmin;
        distancefindAdmin = round(distance*10000);
        dmin10000 = round(dmin*10000);
        mindproperty(1,2) = length(find(distancefindAdmin  == dmin10000)); 
        
%         distancetable = tabulate(distance);
%         dmincountidx = find(distancetable(:,2));
%         mindproperty =  distancetable(dmincountidx,[1,2]);
    end
end
        
