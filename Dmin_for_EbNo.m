function [dmin,Admin,dminforUnitPower] = Dmin_for_EbNo(ridx,mppm,L,EbNoRangeDB,bits,countdmin)
    if size(EbNoRangeDB,1) ==1
        EbNoRangeDB = EbNoRangeDB.';
    end
    [B,Z] = size(mppm);
    signalset = zeros(B*L,Z^L);
%     bits = L*log2(B);
    for ell = 1:L
        signalset(B*(ell-1)+1:B*ell, :) = repmat(kron(mppm,ones(1,Z^(L-ell))) , 1 , Z^(ell-1) );
    end
    if size(ridx,1) ==1
        codebookfull = hadamards(signalset);
        codebook = codebookfull(ridx,:);
    else
        codebook = ridx*signalset;
    end
    [mindproperty,~] = calculateED(codebook,1,countdmin);   % calculate dmin with symbol power per bit = 1;
    dminforUnitPower = mindproperty(1,1);
    if countdmin == 0
        countdmin = size(mindproperty,1);
    end
    EbNoRange = 10.^(EbNoRangeDB./10);
    nNP = EbNoRange.*bits;    % if No = 2, noise variance = 1; then multiply 2 (as mMtc python)
    dmin = mindproperty(:,1) * sqrt(nNP).';  %each coloumn is the dmin for a EbNo
    Admin = repmat(mindproperty(:,2),1,length(EbNoRange));  % each coloumn is the Admin for a EbNo
    dmin = dmin(1:countdmin,:);
    Admin = Admin(1:countdmin,:);
end