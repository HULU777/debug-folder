function mppm = input_(optimalmppm,coding,fitnessmethod,C)
    mppm = ones(1,5)*(-1);  % [Z,B,M,fitness]
    mppm(1:2) = size(optimalmppm);
    mppm(3) = sum(sum(optimalmppm))/mppm(1);
    disp('please be reminded each row is a codeword!!');   
    population = reshape(optimalmppm.',1,[]);
    
    switch coding
        case('Hadamard'); codedpopulation = Hadamardcoding(population);
        case('None'); codedpopulation = population;
    end
    
    switch fitnessmethod
        case 'Q'
            [Ddistribution,mppm(5)] = calculateQ_mind(codedpopulation,mppm(2),mppm(1),C);
            mppm(4) = Ddistribution(1,1);
        case 'square'
            [Ddistribution,mppm(5)] = calculateS_mind(codedpopulation,mppm(2),mppm(1),C,);
            mppm(4) = Ddistribution(1,1);
    end
end