function out = r_con(a, b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%out = 0.25*heaviside(0.5*(a + b - sqrt(a.^2 + b.^2)));
%out = (0.5*(a + b - sqrt(a.^2 + b.^2)));
out = tanh(0.75*(a + b - sqrt(a.^2 + b.^2)));
%out = sin(0.75*(a + b - sqrt(a.^2 + b.^2))./pi);


end

