function [Ddistribution,fitnesssquare] = calculateS_mind(population,B,Z,C,D)  % distance dminparameters
    if nargin<5
        D = 2000;
    end
    a = size(population,1);
    fitnesssquare = zeros(a,1);
    Ddistribution = zeros(a,2);
    for i = 1:a
        f = C;
        codebook = reshape(population(i,:),B,Z);
        [dminproperty,distance] = calculateHD(codebook);
        distance = reshape(distance,1,[]);
        for j = 1:size(distance,2)
            g = distance(1,j);
            if g>D   % if Euclidean, add sqrt()
                g = D;
            end
            f = f - 1/sqrt(g^2+1);
        end
        
%         dminparameters(i,1)  = dmatrix(1,1);
%         dminparameters(i,2) = dmatrix(1,2);
        [dminproperty,~] = sortrows(dminproperty,[1 -2]);
        Ddistribution(i,:) = dminproperty(1,:);
        fitnesssquare(i,1) = f;  %+dminproperty(1,1);
    end
    [dminparametersdisp,di] = sortrows(Ddistribution,[-1 2]);
    disp('maxed_dmin:');disp(dminparametersdisp(1,1));  
    disp('Admin:');disp(dminparametersdisp(1,2));  
    [~,gi] = sort(fitnesssquare);
%     #considering d = 0;
    if gi(1) == di(1)
        disp('square: fitness and dmin ranks the same.');
%         pause;
    end
end