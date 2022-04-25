function cmatrix = HAD16_32()
    hmatrix = hadamard(32);
    cmatrix = hmatrix(15:32,:);
end