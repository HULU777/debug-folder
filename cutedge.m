% the program did not give the best D, 
% because it ends as soon as MPPMd>PPMd;
% delete bigger mind in the matrix not necessarily give better D

clear all; clc;
parameter = [16,16,2,1,8,5];  %[Z,B,W,L,n,lengthM] ;  [42,43,2^6,2];  [2^5,2^5,0,2,400,16,17]
Z = parameter(1);
B = parameter(2);
M = parameter(3);
L = parameter(4);
n= parameter(5);
lengthM = parameter(6);
N = B * L;
EbNo = -10*log10(L*log2(Z));
nrange = 3:15;
nodeset = cell(1,length(nrange));
cmatrix = cell(1,length(nrange));
mppmd = zeros(1,length(nrange));
ppmd = zeros(1,length(nrange));
parfor nn = 1:length(nrange)
    n = nrange(nn);
    [nodeset{nn},cmatrix{nn},mppmd(nn),ppmd(nn)] = vertexsearch(lengthM,Z,B,M,L,N,n,EbNo);
end

function [nMC, cmatrix,MPPMd,PPMd] = vertexsearch(lengthM,Z,B,M,L,N,n,EbNo)
cmatrix = ZCmatrix(lengthM,Z,L,N,n,EbNo,1);
% cmatrix = eye(B);
PPM = MPPMset(B,L,1);
PPMcodebook = cmatrix * PPM;
[mindpropertyPPM,~] = calculateED(PPMcodebook,1,1);
PPMd = mindpropertyPPM(1,1);

MPPM = MPPMset(B,L,M);
MPPMcodebook = cmatrix * MPPM;
[~,~,dmatrix] = calculateED(MPPMcodebook,0,1);
uniqueD = unique(dmatrix);
maxdeletenode = size(MPPM,2) - size(PPM,2);
i = 0; nMC = 0;  lengthD = length(uniqueD);
MPPMd = 0; left = 1; right = lengthD; j = floor((left+right)/2); 
while MPPMd<PPMd && left <= right && (j ~= left || j ~= right)
    ff = setdiff([left right],j);
    if length(ff) == 1 && left +1 == right
        j = ff;
    else
        j = floor((left+right)/2);
    end
    nMCoutput = nMC;
    dmatrix1 = dmatrix<=uniqueD(j);
    dmatrix1 = tril(dmatrix1,-1);
    [row,col] = find(dmatrix1);
    adjacentmatrix = [row';col']';
    nMC=grMinVerCover(adjacentmatrix);
    i = length(nMC);
    if i <= maxdeletenode
        left = j;
        MPPMchosen = MPPM;
        MPPMchosen(:,nMC') = [];
        if size(MPPMchosen,2) > Z
            MPPMchosen1 = MPPMchosen(:,1:Z);
            PPMcodebook = cmatrix * MPPMchosen1;
        else
            PPMcodebook = cmatrix * MPPMchosen;
        end
        [mindpropertyMPPM,~] = calculateED(PPMcodebook,1,1);
        MPPMd = mindpropertyMPPM(1,1);
    else
        right = j-1;
    end
end
end
% MPPMchosen = MPPM;
% MPPMchosen(:,nMCoutput') = [];
% PPMcodebook = cmatrix * MPPMchosen;
% [mindpropertyMPPM,~] = calculateED(PPMcodebook,1,1);
% MPPMd = mindpropertyMPPM(1,1);


