% mppm = [1,2,3;4,5,6;7,8,9;];
% mppm = reshape(mppm,1,[]);
% population = repmat(mppm,2,1);
% population(2,:) = population(2,:) *(-1);
% alpha = 0.3;
% L = 2;
% a = [1    -1     1     1    -1    -1;
%     -1    -1     1    -1     1    -1];
% nM = Hadamardcodings(population,3,3,alpha,L,a);

function [nM,ridx] = Hadamardcoding(population,Z,B,alpha,L,ridx)   % ,table  dmatrix  ,P
% a cloumn of mppm is a codeword
%     P = 10^(P/10);
    pSize = size(population,1);
    population = population.';
    mppm = reshape(population,B,[]);
    signalset = zeros(B*L,pSize*Z^L);
    for p = 1:pSize
        for ell = 1:L
            signalset(B*(ell-1)+1:B*ell, (p-1)*Z^L+1 : p*Z^L) = repmat(kron(mppm(:,(p-1)*Z+1:p*Z),ones(1,B^(L-ell))) , 1 , B^(ell-1) );
        end
    end
    
    % normalize the average power of codeword to P=n, thus, total codebook power is nN. 
    %     Cpower = sum(sum(matrix.^2));   
    %     nM = matrix./ sqrt(Cpower/(P*N));    % normMatrix  n*N
    if ridx == 0
        hadsize = 2^ceil(log2(Z*L));
        hmatrix = hadamard(hadsize);
        n = ceil(alpha*Z*L);
        %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);

        randidxc = randperm(hadsize-1);
        randvalues = values(randidxc);
        columnindex = randvalues(1:Z*L);
        ridx = hmatrix(rowindex,columnindex);
    end
    if size(ridx,1) ==1
        nMfull = hadamards(signalset);
        nM = nMfull(ridx,:);
    else
        nM = cmatrix*signalset;
    end
    
    nM = reshape(nM,[],pSize);
    nM = nM.';
end
        
