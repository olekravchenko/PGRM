%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 2-D cavity problem with
%           Petrov-Galerkin-Rvachev method (PGRM)
%
%               U\theta_X + V\theta_Y = (1/Pr) \grad^2\theta + S(x,y),
%               U\zeta_X  + V\zeta_Y  =        \grad^2\zeta  + Gr \theta_X,
%               \grad^2\psi = -\zeta,
%               U = \psi_Y, V = -\psi_X,          for (x,y) \in \Omega
%               where \Omega domain is defined by R(Rvachev)-functions
%           
%
%              coded by Oleg Kravchenko, BMSTU, 2015.09.29
%              UPD1: 2015.10.11
%              UPD2: 2016.02.24
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] M.A. Basarab. "Solution of stationary two dimensional problems
%     of natural convection in closed domains by R-function method." 
%     Engineering Journal: Science and Innovation (2013): 1-17.  
%     link: http://engjournal.ru/articles/1071/1071.pdf
% [2] M.A. Basarab. "Numerical-analytical method of solving two-dimensional
%     problems of natural convection in a closed cavity."
%     Mathematical Modeling and Computational Methods (2015): 18-35.  
%     link: http://mmcm.bmstu.ru/media/pdf/a40c6baf-f3ea-401c-854d-438f3f2c4517.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Petrov Galerkin Scheme Implementation with the usage of R(Rvachev)-functions method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

addpath 'core'
addpath 'core_rfm'
addpath 'core_basis'

% Initial parameters
L       = 10;
H       = L;
Pr      = 1;
Gr      = 10^3;
Ra      = Pr*Gr;
eps01   = 1e-6;
nx      = 2^6;
ny      = 2^6;
dx      = L / (nx-1);
dy      = H / (ny-1);
tx      = -0.5*L:dx:0.5*L;
ty      = -0.5*H:dy:0.5*H;
[x, y]  = meshgrid(tx, ty);
theta   = 0*x;          % temperature
psi     = 0*x;          % vortex
dzeta   = 0*x;          % flow
u       = 0*x;          % x component of velocity field
v       = 0*x;          % y component of velocity field
% Set of basis functions
nf      = 6;
hf      = 1 / (nf + 1);
% ind_mat = zeros(3,(nf+2)^2);
K       = (nf+2)^2;
Th      = zeros(nx,ny);
Ps      = zeros(nx,ny);
Dz      = zeros(nx,ny);
Th_new  = Th;
Ps_new  = Ps;
Dz_new  = Dz;
ThA     = zeros(K,K);
PsA     = zeros(K,K);
DzA     = zeros(K,K);
ThB     = zeros(1,K);
PsB     = zeros(1,K);
DzB     = zeros(1,K);
BfTh    = zeros(K,nx,ny);
BfPs    = zeros(K,nx,ny);
BfDz    = zeros(K,nx,ny);
ind     = 1:K;

% Omega function
f03 = x + 0.5;
f04 = y + 0.5*H/L;
f05 = 0.5 - x;
f06 = 0.5*H/L - y;
omH = r_con(f03, f04);
omC = r_con(f05, f06);
% Normalized Omega
omH_norm = omH;
omC_norm = omC;
omH_norm(omH<0) = 0;
omC_norm(omC<0) = 0;
fL = omC_norm ./ (omH_norm + omC_norm);
% replace nan's with 0's
fL(isnan(fL)) = 0.5;
% Modified omega
omega = r_con(0.25 - x.^2, 0.25*H/L - y.^2);
omega(omega<0) = 0;
% omega = tanh(30*omega);

% fL = (0.5 + x) .* (0.5 - x) .* (0.25 - y) .* (0.25 + y);
% fL = omega;

% Plot Omega
figure(1),
surf(x, y, fL); view(30, 30);
daspect([1,1,1])
shading interp
light
lighting phong
% material shiny   
% view(0,90)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('BC function $f$','Interpreter','latex')


% Matrix of indices
ind_mat = mat_index(nf);

% compute matrix of basis functions
for k = ind
    i = ind_mat(2,k);
    j = ind_mat(3,k);
    BfTh(k,:,:)   = omega.*bdspln_function(i-1, x, hf, nf).*bdspln_function(j-1, y, hf, nf);
    tmp           = omega.^2.*bdspln_function(i-1, x, hf, nf).*bdspln_function(j-1, y, hf, nf);
    BfPs(k,:,:)   = tmp;
    BfDz(k,:,:)   = del2(tmp, dx, dy);
end

% Matrix of basis functions
% tmp      = bf_mat(ind_mat, K, nx, ny, x, y, hf, nf);
% BfTh(:,:,:) = tmp(:,:,:) .* omega;

res = 1;
%% Computation loop
tic
format long
while res > eps01
    % step 1: Th function
    [fL_x, fL_y]    = gradient(fL,dx,dy);
    [Th_x, Th_y]    = gradient(Th,dx,dy);
    F               = del2(fL,dx,dy)/Pr - (u.*fL_x + v.*fL_y);    
    tmp             = u.*Th_x + v.*Th_y - F;
    % Assembling
    [ThA, ThB]      = pde_assembling(tmp.*omega, BfTh, omega, K, tx, ty, dx, dy);
    % Solution
    ThC             = solve_mat(ThA, ThB);
    % Reconstruct solution
    Th_new          = rec_sol(fL, BfTh, omega, ThC, K);
    
    % step 2: Dz function  
    [Th_x, Th_y]    = gradient(Th_new,dx,dy);
    F               = del2(fL)/Pr - (u.*fL_x + v.*fL_y);    
    tmp             = u.*Th_x + v.*Th_y - F;
    % Assembling
    [DzA, DzB]      = pde_assembling(tmp, BfDz, omega, K, tx, ty, dx, dy);
    % Solution
    PsC             = solve_mat(DzA, DzB);
    % Reconstruct solution
    Dz_new          = rec_sol(fL, BfDz, omega, PsC, K);
    
    % step 3: Psi function
    % Assembling
    [PsA, PsB]      = pde_assembling(-Dz_new, BfPs, omega, K, tx, ty, dx, dy);
    % Solution
    PsC             = solve_mat(PsA, PsB);
    % Reconstruct solution
    Ps_new          = rec_sol(fL, BfPs, omega, PsC, K);
        
    % Update
    Th      = Th_new;
    Dz      = Dz_new;
    Ps      = Ps_new;    
    [u, v]  = gradient(Ps,dx,dy);
    v       = -v;    
    res     = max(max(abs(Th_new - Th)./abs(Th)));
    
    % Output residual
    disp(['res := ' num2str(res)])
end %while

toc;

% plot numerical solution
figure(2),
contourf(x, y, Th)
shading interp
lighting phong
title('Numerical solution $\theta(x,y)$','Interpreter','latex')
axis square

figure(3),
contourf(x, y, Dz)
shading interp
lighting phong
title('Numerical solution $\zeta(x,y)$','Interpreter','latex')
axis square

figure(4),
contourf(x, y, Ps)
shading interp
lighting phong
title('Numerical solution $\psi(x,y)$','Interpreter','latex')
axis square