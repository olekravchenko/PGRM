% test omega_2d
% Sample of the usage of R-operations
% 16/3/16
% Oleg Kravchenko

clc
clear all
close all

addpath '../core_rfm'

a = 1;
tx = linspace(-a, a, 2^7);
ty = tx;

[x, y] = meshgrid(tx, ty);

% params
H = 1;
L = 1;

% used functions
f01 = 0.5^2 - x.^2;
f02 = 0.5^2 - y.^2;

f03 = x + 0.5;
f04 = y + 0.5*H/L;

f05 = 0.5 - x;
f06 = 0.5*H/L - y;

% omegas
om01 = r_diz(f01, f02);
om02 = r_con(f01, f02);
% normalized
om01(om01<0) = 0;
om02(om02<0) = 0;



% fL = omC_norm ./ (omH_norm + omC_norm);

% plot omegas
figure(1),
surf(x, y, om01); view(30, 30);
% daspect([1,1,1])
shading interp
light
lighting phong
% view(0,90)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('BC function $f$','Interpreter','latex')

% plot omegas
figure(2),
surf(x, y, om02); view(30, 30);
% daspect([1,1,1])
shading interp
light
lighting phong
% view(0,90)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('BC function $f$','Interpreter','latex')