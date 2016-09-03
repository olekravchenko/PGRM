function [output1, output2] = bspln_function(m, x, h)

g1    = @(t) t.^3;
g2    = @(t) 1+3*t+3*t.^2-3*t.^3;
d1g1  = @(t) 3*t.^2;
d1g2  = @(t) 3+6*t-9*t.^2;

output1 = 0*x;
output2 = 0*x;

supp1 = (m-2)*h<=x & x<=(m-1)*h;
supp2 = (m-1)*h<=x & x<=(m-0)*h;
supp3 = (m-0)*h<=x & x<=(m+1)*h;
supp4 = (m+1)*h<=x & x<=(m+2)*h;

%function
output1(supp1) = g1(x(supp1) - (m-2)*h)/h^3;
output1(supp2) = g2(x(supp2)/h - (m-1));
output1(supp3) = g2((m+1) - x(supp3)/h);
output1(supp4) = g1((m+2)*h - x(supp4))/h^3;

%1st derivative
output2(supp1) = d1g1(x(supp1) - (m-2)*h)/h^3;
output2(supp2) = d1g2(x(supp2)/h - (m-1))/h;
output2(supp3) = -d1g2((m+1) - x(supp3)/h)/h;
output2(supp4) = -d1g1((m+2)*h - x(supp4))/h^3;

