function [output1, output2] = bdspln_function(m, x, h, nbf)

% output1 = 0*x;
% output2 = 0*x;

%function
switch m
    case 0
        [tmp1, dtmp1] = bspln_function(0,x,h);
        [tmp2, dtmp2] = bspln_function(-1,x,h);
        output1 = tmp1 - 4*tmp2;
        output2 = dtmp1 - 4*dtmp2;
    case 1
        [tmp1, dtmp1] = bspln_function(+1,x,h);
        [tmp2, dtmp2] = bspln_function(-1,x,h);
        output1 = tmp1 - tmp2;
        output2 = dtmp1 - dtmp2;
    case nbf
        [tmp1, dtmp1] = bspln_function(nbf,x,h);
        [tmp2, dtmp2] = bspln_function(nbf+2,x,h);
        output1 = tmp1 - tmp2;
        output2 = dtmp1 - dtmp2;
    case nbf + 1
        [tmp1, dtmp1] = bspln_function(nbf+1,x,h);
        [tmp2, dtmp2] = bspln_function(nbf+2,x,h);
        output1 = tmp1 - 4*tmp2;
        output2 = dtmp1 - 4*dtmp2;
    otherwise
        [output1, output2] = bspln_function(m,x,h);
end %switch

%1st derivative
