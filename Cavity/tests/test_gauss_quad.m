% test Gauss_quad

clc
clear all
close all

addpath '../core'

a = 0; b = pi;
% Gauss quad
i1_1 = Gauss_quad('(sin(x)/x)^2',a,b);
i1_2 = Gauss_quad('(sinc(x))^2',a,b);


f = @(x)(sin(x)./x).^2;
% inline Gauss Krondor
i2_0 = integral(f,0,pi);
i2_1 = quadgk(f,0,pi);
i2_2 = quad(@(x)(sinc(x)).^2,0,pi);
% 
disp([i1_1, i1_2, i2_0, i2_1, i2_2])

%% 2D
format long
a = 3; b = 5;
fun = @(x,y) a*x.^2 + b*y.^2;

[tx, ty] = meshgrid(0:1e-2:5,-5:1e-2:0);

ffun = fun(tx, ty);
q = integral2(ffun,0,5,-5,0,'Method','iterated',...
'AbsTol',0,'RelTol',1e-10);