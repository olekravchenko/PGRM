% test omega_2d
% Sample of the usage of R-operations
% 16/3/16
% Oleg Kravchenko

clc
clear all
close all

addpath '../core_rfm'

a = 8;
tx = linspace(0, a, 2^8);
ty = tx;

[x, y] = meshgrid(tx, ty);

% params
Lx = 8;
Ly = 8;

% used functions
f01 = 1 - (x-1.5).^2 - (y-2).^2;
f01(f01<0) = 0;
f02 = 2 - (x-4).^2 - (y-2).^2;
f02(f02<0) = 0;
f03 = 1 - (x-6.5).^2 - (y-2).^2;
f03(f03<0) = 0;
f04 = 1 - (x-4).^2 - (y-6).^2;
f04(f04<0) = 0;


% omegas
om012       = r_diz(f01, f02);
om0123      = r_diz(om012, f03);
om01234     = r_diz(om0123, f04);

om = tanh(2*om01234);
% fL = omC_norm ./ (omH_norm + omC_norm);

% plot omegas
figure(1),
surf(x, y, om); view(30, 30);
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