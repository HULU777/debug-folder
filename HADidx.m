function cmatrix = HADmatrix(n,N)
    %NO 1st row and column
        hadsize = 2^ceil(log2(N));
        hmatrix = hadamard(hadsize);
        %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);

        randidxc = randperm(hadsize-1);
        randvalues = values(randidxc);
        columnindex = randvalues(1:Z*L);
        cmatrix = hmatrix(rowindex,columnindex);
end