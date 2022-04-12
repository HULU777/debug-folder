function Ddistribution = calculate_Admin(population,B,Z)  % distance dminparameters
    [a,b] = size(population);
    n = b/Z;
    Ddistribution = zeros(a,n*2);
    for i = 1:a
        codebook = reshape(population(i,:),n,Z);
        [dminproperty,~] = calculateED(codebook);
        dminpropertycolumn = reshape(dminproperty.',1,[]);
        Ddistribution(i,1:length(dminpropertycolumn)) = dminpropertycolumn.';%dminproperty(1,:);
    end
    [dminparametersdisp,~] = sortrows(Ddistribution,[-1 2]);  % from good to bad
    disp('maxed_dmin:');disp(dminparametersdisp(1,1));  
    disp('Admin:');disp(dminparametersdisp(1,2));  
end