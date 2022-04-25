% clear;clc;
% L = 2;
% Z = 7;
% hmatrix = hadamard(16);
% cmatrix = hmatrix(5:15,2:15);
% [PPMdmi,PPMAdmin]  = PPMdmins(Z,L,cmatrix);

function [PPMdmin,PPMAdmin] = PPMdmin(Z,L,cmatrix,EbNo)
    ppm = eye(Z);  %BIBD7();
    B = Z;
    signalset = zeros(B*L,Z^L);
    for ell = 1:L
        signalset(B*(ell-1)+1:B*ell, :) = repmat(kron(ppm,ones(1,B^(L-ell))) , 1 , B^(ell-1) );
    end
    codebook = cmatrix*signalset;
    [n,b] = size(codebook);
    EbNo = 10^(EbNo/10);
    nNP = EbNo*log2(b)*2; 
    [mindproperty,~] = calculateED(codebook,nNP);
    PPMdmin = mindproperty(1,1);
    PPMAdmin = mindproperty(1,2);
end