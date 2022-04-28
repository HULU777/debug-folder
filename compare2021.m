countdmin = 0;
EbNoRangedB = -2:0.5:10;
BITS = 12;
load('MPPMSPARC15_512.mat');
MScmatrix = cmatrix;
MSmppm = xbestcl{2};


[MSMPPM_dmin,MSMPPM_Admin,MSdminforUnitPower] = Dmin_for_EbNo(MScmatrix,MSmppm,L,EbNoRangedB,BITS,countdmin);
disp('MSdminforUnitPower:'); disp(MSdminforUnitPower);
MSPe = WEF(MSMPPM_dmin, MSMPPM_Admin,2^(BITS));
semilogy(EbNoRangedB,MSPe); hold on;%  SNRrangedB
xlim([-2,10]); 
ylim([1e-5,1]);

load('SVC_ST(15,1,33,2,1,9).mat');  %svc
svcMPPM_dmin = MPPM_dmin; svcMPPM_Admin = MPPM_Admin; svcdminforUnitPower = dminforUnitPower;
disp('svcdminforUnitPower:'); disp(svcdminforUnitPower);
svcPe = WEF(svcMPPM_dmin, svcMPPM_Admin,2^(BITS));
semilogy(EbNoRangedB,svcPe); hold on;

load('SVC_ST(15,4,9,2,1,5)7303.mat');
stsvcMPPM_dmin = MPPM_dmin; stsvcMPPM_Admin = MPPM_Admin; stsvcdminforUnitPower = dminforUnitPower;
disp('svcdminforUnitPower:'); disp(svcdminforUnitPower);
stsvcPe = WEF(stsvcMPPM_dmin, stsvcMPPM_Admin,2^(BITS));
semilogy(EbNoRangedB,stsvcPe); hold on;

% 9 bits hard to compare
% load('ZC15_512.mat');
% zcMPPM_dmin = MPPM_dmin;  zcMPPM_Admin = MPPM_Admin; zcdminforUnitPower = dminforUnitPower;
% disp('zcdminforUnitPower:'); disp(zcdminforUnitPower);
% zcPe = WEF(zcMPPM_dmin, zcMPPM_Admin,2^(BITS));
% semilogy(EbNoRangedB,zcPe); hold on;

legend('MPPM-SPARC(B=Z=512,L=1,M=7)','SVC(B=33,Z=512,L=1,M=2)','SVC-ST(B=9,Z=32,L=1,M=2,4-QAM)');  % 'ZCsequence-based(B=Z=64,L=2,M=1)'
title('Union bound on 3 different schemes at n=42,k=12.');