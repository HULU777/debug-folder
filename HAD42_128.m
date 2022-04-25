function [cmatrix,rowindex] = HAD42_128(rowindex)
    hadsize = 128;
    hmatrix = hadamard(hadsize);
    n = 42;
    if rowindex ==0
    %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);
    end
    cmatrix = hmatrix(rowindex,1:128);
end