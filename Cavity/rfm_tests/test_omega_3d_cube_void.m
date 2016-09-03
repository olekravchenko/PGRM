% test omega
% Sample of the usage of R-operations
% 12/4/16
% Oleg Kravchenko

addpath '../core_rfm'

clc
% clear all
close all

tx = linspace(-1, 1, 2^8);
ty = tx;
tz = tx;

[x, y, z] = meshgrid(tx, ty, tz);

% used functions
f01 = 0.5 - x;
f02 = x + 0.5;
f03 = 0.5 - y;
f04 = y + 0.5;
f05 = 0.5 - z;
f06 = z + 0.5;

% cutting plane
pA = [0.5 0 0.5];
pB = [0 0.5 0.5];
pC = [0.5 0.5 -0.5];
vAB = pB-pA;
vAC = pC-pA;
vnpln = cross(vAB, vAC);
disp(vnpln)
fpln = vnpln(1)*(x - pA(1)) + vnpln(2)*(y - pA(2)) + vnpln(3)*(z - pA(3));

fsphere = 0.5^2 - (x-0.5).^2 - (y-0.5).^2 - (z-0.5).^2;

% omegas
om01 = r_con(f01, f02);
om02 = r_con(f03, f04);
om03 = r_con(f05, f06);
om04 = r_con(om01, om02);
om05 = r_con(om04, om03);
om06 = r_con(om05, -fsphere);
% om06 = om05;

% normalized omegas
om01_norm = om01;
om02_norm = om02;
om05_norm = om05;
om06_norm = om06;
om01_norm(om01<0) = 0;
om02_norm(om02<0) = 0;
om05_norm(om05<0) = 0;
om06_norm(om06<0) = 0;

% plot omegas
figure(1),
subplot(121)

Sx = -1:.5:1;
Sy = -1:.5:1;
Sz = Sy;
cvals = linspace(-1,1,.5);
% cvals = 1;
contourslice(x,y,z,om06_norm,Sx,Sy,Sz,cvals)
axis([-1,1,-1,1,-1,1])
daspect([1,1,1])
view(150, 30)
% campos([0,-20,7])
box on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
% title('Cube with cutted edge defined by $\omega$ function','Interpreter','latex')
title('Conotour plot of $\omega$ function','Interpreter','latex')

% figure(2),

subplot(122)
p = patch(isosurface(x,y,z,om06_norm,1e-6));
isonormals(x,y,z,om06_norm,p)
purple = [1.0 0.5 1.0];
p.FaceColor = purple;
p.EdgeColor = 'none';
axis([-1,1,-1,1,-1,1])
daspect([1,1,1])
view(150, 30);
camlight left
lighting gouraud
box on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
% title('Cube with cutted edge defined by $\omega$ function','Interpreter','latex')
title('Solid plot of $\omega$ function','Interpreter','latex')