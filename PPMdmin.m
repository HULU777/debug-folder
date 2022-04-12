% clear;clc;
% L = 2;
% Z = 7;
% hmatrix = hadamard(16);
% cmatrix = hmatrix(5:15,2:15);
% [PPMdmi,PPMAdmin]  = PPMdmins(Z,L,cmatrix);

function [PPMdmin,PPMAdmin] = PPMdmin(Z,L,cmatrix)
    ppm = eye(Z);
    B = Z;
    signalset = zeros(B*L,Z^L);
    for ell = 1:L
        signalset(B*(ell-1)+1:B*ell, :) = repmat(kron(ppm,ones(1,B^(L-ell))) , 1 , B^(ell-1) );
    end
    codebook = cmatrix*signalset;
    [mindproperty,~] = calculateE2D(codebook);
    PPMdmin = mindproperty(1,1);
    PPMAdmin = mindproperty(1,2);
end