function cmatrix = HADCidx(n,N)
    hmatrix = hadamard(N);
    % real part
    values = 2:N;
    randidxr = randperm(N-1);
    randvalues = values(randidxr);
    rowindex = randvalues(1:n);
    % imaginary part
    randidxrC = randperm(N-1);
    randvaluesC = values(randidxrC);
    rowindexC = randvaluesC(1:n);
    % MATRIX
    cmatrix = hmatrix(rowindex,:)+ hmatrix(rowindexC,:)*1i;
end