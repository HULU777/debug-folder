close all; clear all;
lengthM = 7;
k= 0:lengthM-1;
F = zeros(lengthM,lengthM*(lengthM-1));
for u = 1: lengthM-1
    for e=0:lengthM-1
        idx = (u-1)*lengthM+e+1;
        shiftk = mod(k+e,lengthM);
        F(:, idx) = exp(-1i.*pi./lengthM.*u.*shiftk.*(shiftk+mod(lengthM,2)));
    end
end

correlation = 100 * ones(size(F,2));
for i = 1:size(F,2)
    for j = i + 1:size(F,2)
        correlation(i,j) = F(:,i)'*F(:,j);
    end
end
aa1 = real(correlation);
aa1 = reshape(aa1,1,[]);
RE = sort(aa1);

sizeF2 = size(F,2);
Fmppm = zeros(lengthM,nchoosek(sizeF2,2));
mppmidx = nchoosek(1:sizeF2,2);
k = 1;
while k<=nchoosek(sizeF2,2)
    Fmppm(:,k) = F(:,mppmidx(k,1)) + F(:,mppmidx(k,2));
    k=k+1;
end
mppmcorrelation = 100 * ones(size(Fmppm,2));
for i = 1:size(Fmppm,2)
    for j = i + 1:size(Fmppm,2)
        mppmcorrelation(i,j) = Fmppm(:,i)'*Fmppm(:,j);
    end
end

aa = real(mppmcorrelation);
aa = reshape(aa,1,[]);
mppmRE = sort(aa);