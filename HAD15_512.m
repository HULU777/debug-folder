function [cmatrix,rowindex] = HAD15_512(rowindex)
    hadsize = 512;
    hmatrix = hadamard(hadsize);
    n = 15;
    if rowindex ==0
    %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);
    end
    cmatrix = hmatrix(rowindex,1:hadsize);
end