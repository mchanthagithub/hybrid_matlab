function [root] = secondRootOfQuadratic(a,b,c,dscr_sqrt)
if(b > 0)
    root = 2*c / (-b - dscr_sqrt);
else
    root = (-b + dscr_sqrt) / (2.0 * a);
end


end

