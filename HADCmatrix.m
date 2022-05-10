function cmatrix = HADCmatrix(n,N)
        hadsize = 2^ceil(log2(max(n+1,N)));
        hmatrix = hadamard(hadsize);
% real part
        %NO 1st row
        values = 2:hadsize;
        randidxr = randperm(hadsize-1);
        randvaluesr = values(randidxr);
        rowindex = randvaluesr(1:n);

        randidxc = randperm(hadsize);
        columnindex = randidxc(1:N);

% imaginary part
        randidxi1 = randperm(hadsize-1);
        randvaluesri1 = values(randidxi1);
        rowindexi1 = randvaluesri1(1:n);

        randidxci1 = randperm(hadsize);
        columnindexi1 = randidxci1(1:N);
    % MATRIX
    cmatrix = hmatrix(rowindex,columnindex)+ hmatrix(rowindexi1,columnindexi1)*1i;
end