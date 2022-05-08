function ridx = ZCmatrix(lengthM,Z,L,N,n,EbNo,samples)
k= 0:lengthM-1;
F = zeros(lengthM,lengthM*(lengthM-1));
for u = 1: lengthM-1
    for e=0:lengthM-1
        idx = (u-1)*lengthM+e+1;
        shiftk = mod(k+e,lengthM);
        F(:, idx) = exp(-1i.*pi./lengthM.*u.*shiftk.*(shiftk+mod(lengthM,2)));
    end
end

i = 1; dminforUnitPowerbest = 0;
while i <=samples
    cidx = randperm(lengthM*(lengthM-1));
    selectcidx = cidx(1:N);
    deleteridx = ceil(rand(lengthM-n,1)*lengthM);
    cmatrix = F(:,selectcidx);
    cmatrix(deleteridx,:) = [];

    [dminforUnitPower,~] = PPMdmin(Z,L,cmatrix,EbNo);
    if dminforUnitPower > dminforUnitPowerbest
        selectcidxbest = selectcidx;
        deleteridxbest = deleteridx;
        dminforUnitPowerbest = dminforUnitPower;
%         break;
    end
    i = i+1;
end
ridx = F(:,selectcidxbest);
ridx(deleteridxbest,:) = [];
end