function rowindex = HADidx(n,N)
    %NO 1st row and column
        values = 2:N;
        randidxr = randperm(N-1);
        randvalues = values(randidxr);
        rowindex = randvalues(1:n);
end