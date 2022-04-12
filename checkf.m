function [maxf, bestIndexf,check] = checkf(fitness,optimalf)
    check = 0;
    [maxf, bestIndexf] = max(fitness);
    if (maxf>=optimalf)
        check = 1;
        disp('fitness achieved.');
    end
end