%% test EbNo/SNR
clear all; close all;
EbNoRangedB = 0:10;
cmatrix = qammod(0:15, 16);
countdmin = 0;
L=1;
[MPPM_dmin,MPPM_Admin] = Dmin_for_EbNo(cmatrix,eye(16),L,EbNoRangedB,countdmin);
Pe = WEF(MPPM_dmin, MPPM_Admin,16);
semilogy(EbNoRangedB,Pe);
xlim([0,10]); 
ylim([1e-4,1]);