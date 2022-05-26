

















% left = 1; right = 10;
% while left <= right && (j ~= left || j ~= right)
%     ff = setdiff([left right],j);
%     if length(ff) == 1 && left +1 == right
%         j = ff;
%     else
%         j = floor((left+right)/2);
%     end
%    
%     if j <= 1
%         left = j;
%     else
%         right = j-1;
%     end
% end