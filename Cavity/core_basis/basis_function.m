function [output] = basis_function(k, x, y)

l = floor(.5*(sqrt(1+8*(k-1))-1));
j = k - .5*l*(l+1) - 1;
i = l - j;

output = x.^i .* y.^j;