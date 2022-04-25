function Ddistribution = calculate_Admin(population,B,Z,EbNo)  % distance dminparameters
    [a,b] = size(population);
    n = b/Z;   % Z codewords in each codebook
    Ddistribution = zeros(a,n*2);
    EbNo = 10^(EbNo/10);
    nNP = EbNo*log2(Z)*2;   % Power of the whole codebook sampling theorem
    for i = 1:a
        codebook = reshape(population(i,:),n,Z);
        [dminproperty,~] = calculateED(codebook,nNP);
        dminpropertycolumn = reshape(dminproperty.',1,[]);
        Ddistribution(i,1:length(dminpropertycolumn)) = dminpropertycolumn.';%dminproperty(1,:);
    end
    [dminparametersdisp,~] = sortrows(Ddistribution,[-1 2]);  % from good to bad
    disp('maxed_dmin:');disp(dminparametersdisp(1,1));  
    disp('Admin:');disp(dminparametersdisp(1,2));  
end