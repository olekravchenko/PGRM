% test omega

clc
clear all
close all

addpath '../core_rfm'

a = 0.5;
tx = linspace(-a, a, 512);
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

omH = r_diz(f03, f04);
omC = r_diz(f05, f06);

% normalized omegas
om01_norm = om01;
om02_norm = om02;
om01_norm(om01<0) = 0;
om02_norm(om02<0) = 0;

omH_norm = omH;
omC_norm = omC;
omH_norm(omH<0) = 0;
omC_norm(omC<0) = 0;

% fL = omC_norm ./ (omH_norm + omC_norm);

% plot omegas
figure(1),
surf(x, y, omH_norm); view(30, 30);
daspect([1,1,1])
shading interp
light
lighting phong
% view(0,90)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('BC function $f$','Interpreter','latex')