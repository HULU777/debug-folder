close all; clear all;
% X axis: SNR 
parameter =  [42,43,2^6,2];  %[n,lengthM,B,L]   [16,17,2^8,1]  [42,43,2^6,2]
n=parameter(1); % 64
lengthM=parameter(2); %67
B = parameter(3);  % 2^12
Z = B;
L = parameter(4);
N = B*L;
bits = L*log2(B);
EbNo = 1;
countdmin = 1;
k= 0:lengthM-1;
F = zeros(lengthM,lengthM*(lengthM-1));
for u = 1: lengthM-1
    for e=0:lengthM-1
        idx = (u-1)*lengthM+e+1;
        shiftk = mod(k+e,lengthM);
        F(:, idx) = exp(-1i.*pi./lengthM.*u.*shiftk.*(shiftk+mod(lengthM,2)));
    end
end

SNRrangedB = -6:1:4;
EbNoRangeSNRdB = SNRrangedB - 10* log10( bits);

EbNoRangedB = -2:0.2:10;
i = 1; dminforUnitPowerbest = 0.9070;
while i <=500
    cidx = randperm(lengthM*(lengthM-1));
    selectcidx = cidx(1:N);
    deleteridx = ceil(rand(lengthM-n)*lengthM);
    cmatrix = F(:,selectcidx);
    cmatrix(deleteridx,:) = [];

    [MPPM_dmin,MPPM_Admin,dminforUnitPower] = Dmin_for_EbNo(cmatrix,eye(B),L,EbNoRangedB,countdmin); % EbNoRangedB  EbNoRangeSNRdB
    if dminforUnitPower > dminforUnitPowerbest
        MPPM_dminbest = MPPM_dmin;
        MPPM_Adminbest = MPPM_Admin;
        cidxbest = cidx;
        deleteridxbest = deleteridx;
        dminforUnitPowerbest = dminforUnitPower;
        break;
    end
    i = i+1;
end
disp('dminforUnitPower:'); disp(dminforUnitPowerbest);
Pe = WEF(MPPM_dminbest, MPPM_Adminbest,Z^L);
semilogy(EbNoRangedB,Pe); %  SNRrangedB
xlim([-2,10]); 
ylim([1e-5,1]);