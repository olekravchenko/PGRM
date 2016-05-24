function [a, b, c, d, x, y] = RegionIC(id, varargin)

nVarargs = length(varargin);

if nVarargs == 0
    nx = 2^16;
else
    disp('Incorrect Input');
    a = 0; b = 0; c = 0; d = 0;  
    return
end % if

if id == 1
    a = -1; b = -a;
    c = -1; d = -c;
else
    disp('Incorrect Input');
    a = 0; b = 0; c = 0; d = 0;  
    return
end % if