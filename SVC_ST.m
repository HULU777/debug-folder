% close all; clear all;
%% parameter and coding matrix
% [15,4,9,2,1,5]
% n=16; k=8,10bits   [16,4,7,2,1,4]  [16,4,12,2,1,6]
parameter =  [16,4,7,2,1,4] ;  %[n,QAM,B,M,L,b0]    [15,1,33,2,1,9]    [15,1,9,2,1,5]
n=parameter(1); % 64
QAM=parameter(2); %67
B = parameter(3);  % 2^12
M = parameter(4);
L = parameter(5);
b0 = parameter(6);
Z = 2^b0;
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
idx = randperm(nchoosek(B,M));  nonzeroposition = nonzeroposition(idx(1:2^b0),:);
G = zeros(M,QAM); sparset = zeros(B,2^b0,2);
for column = 1:M
    if QAM ==1
        G(1,:) = 1; G(2,:) = 1j;
    else
        G(column,:) = qammod((0:QAM-1).',QAM) * exp(1i*(pi*(column-1)/(2*M)));
    end
    for i = 1: 2^b0
        nzposition = nonzeroposition(i,column);
        sparset(nzposition,i,column) = 1;
    end
end
sparsetQAM1 = kron(sparset(:,:,1),repmat(G(1,:),1,QAM));
sparsetQAM2= kron(sparset(:,:,2),kron(G(2,:),ones(1,QAM)));
sparsetQAM = sparsetQAM1 + sparsetQAM2;

bits = (b0+M*log2(QAM))*L;
countdmin = 0;

SNRrangedB = 0:8;
EbNoRangeSNRdB = SNRrangedB - 10*log10( bits);

EbNoRangedB = -2:0.2:12;

k = 0; dminforUnitPowerbest = 0;
while k < 10
    miu = -1;
    cmatrix = rand(n,N)>0.5;
    cmatrix = -1*cmatrix;
    cmatrix(cmatrix==0) = 1;
    hmatrix = hadamard(32);
%     cmatrix = hmatrix(7:22,2:13);  % 4:23,2:10
    cmatrix = hmatrix(7:22,7:13);
    for i = 1:n-1
        for j = i+1:n
            miu = max(miu,cmatrix(i,:) * cmatrix(j,:).'/n);
        end
    end
    disp(miu);
    [MPPM_dmin,MPPM_Admin,dminforUnitPower] = Dmin_for_EbNo(cmatrix,sparsetQAM,L,EbNoRangedB,bits,countdmin); % EbNoRangedB  EbNoRangeSNRdB
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
Pe = WEF(MPPM_dmin, MPPM_Admin,2^(b0+M*log2(QAM)));
semilogy(EbNoRangedB,Pe); hold on;%  SNRrangedB
xlim([-2,12]); 
ylim([1e-5,1]);