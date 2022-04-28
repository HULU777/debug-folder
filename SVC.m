close all; clear all;
%% parameter and coding matrix
parameter =  [15,92,2,1,12];  %[n,B,M,L,b] [42,33,2,1,9] [42,92,2,1,12]
n=parameter(1); % 64
B = parameter(2);  % 2^12
Z = B;
M = parameter(3);
L = parameter(4);
b = parameter(5);
N = B*L;

%% sparse vector set
% nonzeroposition = nchoosek(1:B,M);
% idx = randperm(nchoosek(B,M));  idx = nonzeroposition(idx(1:2^b0));
% sparsetQAM = zeros(B,2^b0); G = zeros(M,QAM);
% for column = 1:M
%     sparset = zeros(B,2^b0);
%     G(column,:) = qammod((0:3).',4) * exp(1i*(pi*(column-1)/(2*M)));
%     for i = 1: 2^b0
%         nzposition = nonzeroposition(i,column);
%         sparset(nzposition,i) = 1;
%     end
%     sparsetQAM= sparsetQAM + kron(sparset,G(column,:));
% end
% sparsetQAM = sparsetQAM(:,:,column) + sparsetQAM(:,:,column);
% sparset = repmat(sparset,1,QAM^M);


nonzeroposition = nchoosek(1:B,M);
idx = randperm(nchoosek(B,M));  nonzeroposition = nonzeroposition(idx(1:2^b),:);
sparset = zeros(B,2^b);
for i = 1: 2^b
    nzposition = nonzeroposition(i,:);
    sparset(nzposition,i) =[ 1, 1j];
end

bits = b*L;
countdmin = 0;

SNRrangedB = 0:8;
EbNoRangeSNRdB = SNRrangedB - 10*log10( bits);

EbNoRangedB = -2:10;

k = 0; dminforUnitPowerbest = 0.6220;
while k < 100
    miu = -1;
    cmatrix = rand(n,N)>0.5;
    cmatrix = -1*cmatrix;
    cmatrix(cmatrix==0) = 1;
%     hmatrix = hadamard(64);
%     cmatrix = hmatrix(7:21,2:34);  % 4:23,2:10
    for i = 1:n-1
        for j = i+1:n
            miu = max(miu,cmatrix(i,:) * cmatrix(j,:).'/n);
        end
    end
    disp(miu);
    [MPPM_dmin,MPPM_Admin,dminforUnitPower] = Dmin_for_EbNo(cmatrix,sparset,L,EbNoRangedB,bits,countdmin); % EbNoRangedB  EbNoRangeSNRdB
    if dminforUnitPower > dminforUnitPowerbest
        MPPM_dminbest = MPPM_dmin;
        MPPM_Adminbest = MPPM_Admin;
        cmatrixbest = cmatrix;
        dminforUnitPowerbest = dminforUnitPower;
        break;
    end
    k = k+1;
end

disp('dminforUnitPowerbest:'); disp(dminforUnitPowerbest);
Pe = WEF(MPPM_dminbest, MPPM_Adminbest,2^(bits));
semilogy(EbNoRangedB,Pe); %  SNRrangedB
xlim([-2,10]); 
ylim([1e-5,1]);