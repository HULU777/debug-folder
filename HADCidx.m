function cmatrix = HADCidx(n,N)
        hadsize = 2^ceil(log2(N));
        hmatrix = hadamard(hadsize);
% real part
        %NO 1st row and column
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvaluesr = values(randidxr);
        rowindex = randvaluesr(1:n);

        randidxc = randperm(hadsize-1);
        randvaluesc = values(randidxc);
        columnindex = randvaluesc(1:Z*L);

% imaginary part
        randidxi1 = randperm(hadsize-1);
        randvaluesri1 = values(randidxi1);
        rowindexi1 = randvaluesri1(1:n);

        randidxci1 = randperm(hadsize-1);
        randvaluesci1 = values(randidxci1);
        columnindexi1 = randvaluesci1(1:Z*L);
    % MATRIX
    cmatrix = hmatrix(rowindex,columnindex)+ hmatrix(rowindexi1,columnindexi1)*1i;
end