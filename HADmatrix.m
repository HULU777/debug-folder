function cmatrix = HADmatrix(n,N)
        hadsize = 2^ceil(log2(max(n+1,N)));
        hmatrix = hadamard(hadsize);
        %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);

        randidxc = randperm(hadsize);
        columnindex = randidxc(1:N);
        cmatrix = hmatrix(rowindex,columnindex);
end