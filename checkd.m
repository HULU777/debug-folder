function [mind, Admin, bestIndexd,check] = checkd(Ddistribution,optimald)
    check = 0;
    [mindmatrix, sortedIndexd] = sortrows(Ddistribution,[-1 2]); % from good to bad
    mind = mindmatrix(1,1);
    Admin = mindmatrix(1,2);
    if (mind>=optimald)
        check = 1;
        disp('DMIN achieved.');
    end
    bestIndexd = sortedIndexd(1);
end