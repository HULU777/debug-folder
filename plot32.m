%% given Dmin and Admin at EbNo=1, plot for different EbNos 
function plot32(bits, EbNoRange,Z,L)
    if size(EbNoRange,1)==1
        EbNoRange = EbNoRange.';
    end
%     load('DandA(B=16,L=2)');
%     My.L = 2; My.Z = 16;
%     count = My.Z^(My.L);
%     count = count*(count-1);
    EbNo = 10.^(EbNoRange./10);
    nNP = EbNo.*bits.*2; 
    for i = 1:length(D)
        dmin = D(i) * sqrt(nNP).';  %each coloumn is the dmin for a EbNo
        Admin = repmat(A,1,length(EbNo));  % each coloumn is the Admin for a EbNo
        Pe = WEF(dmin,Admin,Z^L);
        semilogy(EbNodb,Pe);hold on;
    end
    title('Z=B=16,L=2,n=32,R=0.25');
end