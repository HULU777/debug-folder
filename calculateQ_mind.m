function [Ddistribution,fitnessQ] = calculateQ_mind(population,B,Z,C,Qdominant)  % distance dminparameters
    a = size(population,1);
    fitnessQ = zeros(a,1);
    Ddistribution = zeros(a,2);
    for i = 1:a
        fq = C;
        codebook = reshape(population(i,:),B,Z);
        [dminproperty,distance] = calculateHD(codebook);
        distance = reshape(distance,1,[]);
        for j = 1:size(distance,2)
            g = distance(1,j);
%             fq = fq - 1/sqrt(g^2+1);
            fq = fq - qfunc(g/Qdominant);%6.6667);  %3.8911
        end
        
%         dminparameters(i,1)  = dmatrix(1,1);
%         dminparameters(i,2) = dmatrix(1,2);
        [dminproperty,~] = sortrows(dminproperty,[1 -2]);
        Ddistribution(i,:) = dminproperty(1,:);
        fitnessQ(i,1) = fq;
    end
    [dminparametersdisp,di] = sortrows(Ddistribution,[-1 2]);
    disp('maxed_dmin:');disp(dminparametersdisp(1,1));  
    disp('Admin:');disp(dminparametersdisp(1,2));  
    [~,giq] = sort(fitnessQ);
    if giq(1) == di(1)
        disp('Q: fitness and dmin ranks the same.');
%         pause;
    end
end