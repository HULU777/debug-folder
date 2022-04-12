% mppm = BIBD7();
% alpha = 0.3;
% L = 4;
% nM = Hadmardcodings(mppm,alpha,L);

function nM = Hadmardcoding(mppm,alpha,L)   % ,table  dmatrix  ,P
% a cloumn of mppm is a codeword
%     P = 10^(P/10);
    [B,N] = size(mppm);
    signalset = zeros(B*L,N^L);
    for ell = 1:L
        signalset(B*(ell1)+1:B*ell, :) = repmat(kron(mppm,ones(1,B^(L-ell))) , 1 , B^(ell-1) );
    end
    
    % normalize the average power of codeword to P=n, thus, total codebook power is nN. 
    %     Cpower = sum(sum(matrix.^2));   
    %     nM = matrix./ sqrt(Cpower/(P*N));    % normMatrix  n*N

    hadsize = 2^ceil(log2(N*L));
    hmatrix = hadamard(hadsize);
    n = ceil(alpha*N*L);
    %NO 1st row and column
    values = 2:hadsize;
    randidxr = randperm(hadsize-1);
    randvalues = values(randidxr);
    rowindex = randvalues(1:n);

    randidxc = randperm(hadsize-1);
    randvalues = values(randidxc);
    columnindex = randvalues(1:N*L);
    cmatrix = hmatrix(rowindex,columnindex);
    nM = cmatrix*signalset;
end
        
