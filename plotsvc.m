close all; clear all;
parameter1 = [15,4,9,2,1,5];
% parameter2 = [15,16,9,2,2,5];
parameter3 = [15,33,2,1,9];
% parameter4 = [15,129,2,1,13];
EbNoRangedB = -2:14;
Pe1 =  plotsvcst(parameter1,EbNoRangedB);
% Pe2 =  plotsvcst(parameter2,EbNoRangedB);
Pe3 =  plotsvcs(parameter3,EbNoRangedB);
% Pe4 =  plotsvcs(parameter4,EbNoRangedB);
semilogy(EbNoRangedB,Pe1,'*-','Color','[1 0.5 0]'); hold on;%  SNRrangedB
% semilogy(EbNoRangedB,Pe2,'o-','Color','[1 0.5 0]'); hold on;
semilogy(EbNoRangedB,Pe3,'b*-'); hold on;
% semilogy(EbNoRangedB,Pe4,'bo-'); hold on;
xlim([-2,14]); 
ylim([1e-5,1]);
legend("SVC-ST,b_t=9,bn=2","SVC,b_s=9")  %"SVC-ST,b_t=13,bn=4", ,"SVC,b_s=13"

function Pe =  plotsvcs(parameter,EbNoRangedB)  %[n,B,M,L,b] [42,33,2,1,9]   [15,92,2,1,12]
n=parameter(1); % 64
B = parameter(2);  % 2^12
Z = B;
M = parameter(3);
L = parameter(4);
b = parameter(5);
N = B*L;


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


k = 0; dminforUnitPowerbest = 0;
while k < 500
    miu = -1;
    cmatrix = rand(n,N)>0.5;
    cmatrix = -1*cmatrix;
    cmatrix(cmatrix==0) = 1;
%     hmatrix = hadamard(64);
%     cmatrix = hmatrix(7:21,10:42);  % 4:23,2:10
    for i = 1:n-1
        for j = i+1:n
            miu = max(miu,cmatrix(i,:) * cmatrix(j,:).'/n);
        end
    end
%     disp(miu);
    [MPPM_dmin,MPPM_Admin,dminforUnitPower] = Dmin_for_EbNo(cmatrix,sparset,L,EbNoRangedB,bits,countdmin); % EbNoRangedB  EbNoRangeSNRdB
    if dminforUnitPower > dminforUnitPowerbest
        MPPM_dminbest = MPPM_dmin;
        MPPM_Adminbest = MPPM_Admin;
        cmatrixbest = cmatrix;
        dminforUnitPowerbest = dminforUnitPower;
%         break;
    end
    k = k+1;
end
disp('dminforUnitPowerbest:'); disp(dminforUnitPowerbest);
Pe = WEF(MPPM_dminbest, MPPM_Adminbest,2^(bits));
end

function Pe = plotsvcst(parameter,EbNoRangedB)
n=parameter(1); % 64
QAM=parameter(2); %67
B = parameter(3);  % 2^12
M = parameter(4);
L = parameter(5);
b0 = parameter(6);
Z = 2^b0;
N = B*L;

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

k = 0; dminforUnitPowerbest = 0;
while k < 1
    miu = -1;
    cmatrix = rand(n,N)>0.5;
    cmatrix = -1*cmatrix;
    cmatrix(cmatrix==0) = 1;
    hmatrix = hadamard(32);
    cmatrix = hmatrix(7:21,2:10);  % 4:23,2:10
    for i = 1:n-1
        for j = i+1:n
            miu = max(miu,cmatrix(i,:) * cmatrix(j,:).'/n);
        end
    end
%     disp(miu);
    [MPPM_dmin,MPPM_Admin,dminforUnitPower] = Dmin_for_EbNo(cmatrix,sparsetQAM,L,EbNoRangedB,bits,countdmin); % EbNoRangedB  EbNoRangeSNRdB
    if dminforUnitPower > dminforUnitPowerbest
        MPPM_dminbest = MPPM_dmin;
        MPPM_Adminbest = MPPM_Admin;
        cmatrixbest = cmatrix;
        dminforUnitPowerbest = dminforUnitPower;
%         break;
    end
    k = k+1;
end
disp('dminforUnitPowerbest:'); disp(dminforUnitPowerbest);
Pe = WEF(MPPM_dminbest, MPPM_Adminbest,2^(b0+M*log2(QAM)));
end
