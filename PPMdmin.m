% clear;clc;
% L = 2;
% Z = 7;
% hmatrix = hadamard(16);
% cmatrix = hmatrix(5:15,2:15);
% [PPMdmi,PPMAdmin]  = PPMdmins(Z,L,cmatrix);

function [PPMdmin,PPMAdmin] = PPMdmin(Z,L,ridx,EbNo)
    ppm = eye(Z);  %BIBD7();
    B = Z;
    signalset = zeros(B*L,Z^L);
    for ell = 1:L
        signalset(B*(ell-1)+1:B*ell, :) = repmat(kron(ppm,ones(1,B^(L-ell))) , 1 , B^(ell-1) );
    end
    if size(ridx,1) ==1
        codebookfull = hadamards(signalset);
        codebook = codebookfull(ridx,:);
    else
        codebook = cmatrix*signalset;
    end
    [n,b] = size(codebook);
    EbNo = 10^(EbNo/10);
    nNP = EbNo*log2(b); 
    [mindproperty,~] = calculateED(codebook,nNP,1);
    PPMdmin = mindproperty(1,1);
    PPMAdmin = mindproperty(1,2);
end