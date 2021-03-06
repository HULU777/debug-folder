function [t, mind, dMatrixMPPM, MPPM] = MIMO(B,L,M,Z,cmatrix,EbNoDB)
    tic
    i = 0;
    MPPM = MPPMset(B,L,M);
    while size(MPPM,2) > Z^L
        codebook = cmatrix * MPPM;
        z = size(codebook,2);
        EbNodb = EbNoDB *log10(log2(z)) /log10(log2(Z^L));
        EbNo = 10^(EbNodb/10);
        nNPinput =  EbNo*log2(z);
        [~,dMatrix] = calculateED(codebook,0,2);
        while i == 0
            dMatrixMPPM = dMatrix;
            i =1;
        end
        dMatrixsort =  sort(dMatrix,2,'descend');
        dMatrixsort = dMatrixsort(:,1:Z-1);
        [A,  minI] = sortrows(dMatrixsort, flip(1:Z-1));
        MPPM(:,minI(1)) = [];
    end
    codebook = cmatrix * MPPM;
    z = size(codebook,2);
    EbNodb = EbNoDB *log10(log2(z)) /log10(log2(Z^L));
    EbNo = 10^(EbNodb/10);
    nNPinput = EbNo*log2(z);
    [~,dMatrix] = calculateED(codebook,0,2);
    while i == 0
        dMatrixMPPM = dMatrix;
        i =1;
    end
    mind = min(min(dMatrix));
    t = toc;
end


%         [x,y]=find(dMatrix==minimum);
%         if y(1)>=x(1)
%             y(1) = y(1)+1;
%         end
%         dx = sort(dMatrix(x(1),:));
%         dy = sort(dMatrix(y(1),:));
%         if dx(2) > dy(2)
%             MPPM(:,y(1)) = [];
%         else
%             MPPM(:,x(1)) = [];
%         end