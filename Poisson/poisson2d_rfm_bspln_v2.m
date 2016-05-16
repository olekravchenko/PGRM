%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 2-D Poisson equation with
%           Petrov-Galerkin-Rvachev method (PGRM)
%
%               \nabla^2 u(x,y) = f(x,y),
%               u(x,y)|_{\partial\Omega}  = 0,              (Dirichlet problem)
%               where \Omega domain is defined by R(Rvachev)-functions
%           
%           Here u^N(x,y) = \sum_k^N \psi_k(x,y), \psi_k = B_i(x)B_j(y), 
%               k is multiindex, i=i(k), j=j(k).
%               B_i(x)B_j(y) is a tensor product of 1D cubic B-splines,
%               from [4] .
%
%              coded by Oleg Kravchenko, BMSTU, 2015.07.02
%              UPD1: 2015.10.01
%              UPD2: 2016.04.02
%              UPD3: 2016.05.16
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] V.F. Kravchenko, M.A. Basarab. "Boolean Algebra and Approximation
%     Methods in Boundary-Value Problems of Electrodynamic. Fizmatlit,
%     Moscow, 2004. Pages, 308. [in Russian]
%     link: http://gen.lib.rus.ec/book/index.php?md5=1C3E1738CC9BC8C0B90FB63F5968EED0
% [3] M.A. Basarab. "Numerical-analytical method of solving two-dimensional
%     problems of natural convection in a closed cavity."
%     Mathematical Modeling and Computational Methods (2015): 18-35.  
%     link: http://mmcm.bmstu.ru/media/pdf/a40c6baf-f3ea-401c-854d-438f3f2c4517.pdf
% [2] M.A. Basarab. "Solution of stationary two dimensional problems
%     of natural convection in closed domains by R-function method." 
%     Engineering Journal: Science and Innovation (2013): 1-17.  
%     link: http://engjournal.ru/articles/1071/1071.pdf
% [3] M.A. Basarab. "Numerical-analytical method of solving two-dimensional
%     problems of natural convection in a closed cavity."
%     Mathematical Modeling and Computational Methods (2015): 18-35.  
%     link: http://mmcm.bmstu.ru/media/pdf/a40c6baf-f3ea-401c-854d-438f3f2c4517.pdf
% [4] S. Beck. Cubic spline finite element method for solving Poisson's equation
%     on a square. M.Sc. Thesis. Colorado School of Mines. 2008/
%     link: https://dspace.library.colostate.edu/bitstream/handle/11124/353/Beck_mines_0052N_10391.pdf?sequence=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Petrov Galerkin Scheme Implementation with the usage of R(Rvachev)-functions method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Pathes
addpath 'core' 'core_rfm'
%% Initial parameters
a = 0; b = 5;
c = 0; d = 1;
% Space grid
nx = 2^5;   dx = (b-a) / (nx-1);
ny = 2^5;   dy = (d-c) / (ny-1);
% 
tx = linspace(a,b,nx);
ty = linspace(c,d,ny);
[x, y] = meshgrid(tx, ty);
% 
nf  = 3;
hf  = 1 / (nf + 1);
K   = (nf+2)^2;
ind = 1:K;
iC  = 1:min([nx ny]);
% 
id_omega    = 3;                         % omega function
id_bc       = 1;                         % BC, r.h.s functions


%% Load BC and r.h.s. function
[bc_mat, rhs_mat, gf_mat] = init_mat(id_bc, x, y);
%% Take omega function
omega = omega_mat(id_omega, x, y);
% Normalization
omega(omega<0)      = 0;
gf_mat(gf_mat<0)    = 0;
% Modification
omega   = tanh(100*omega);


figure(1),
surf(x, y, omega)
shading interp
%lighting phong
title('\omega_{{w}}(x,y)=tanh(300 \omega(x,y))','Interpreter','tex')
axis equal
%camlight('left');

% Show basis
figure(2),
show_basis(omega, nf, hf, x, y)

tic
%% Matrix of indices
ind_mat = mat_index(nf);
%% Matrix of basis functions
BF      = bf_mat(ind_mat, K, nx, ny, x, y, hf, nf);
%% Assembling
[A, B]  = pde_assembling(rhs_mat, BF, omega, K, tx, ty, dx, dy);
%% Solution
C       = solve_mat(A, B);
%% Reconstruct solution
u_appr  = rec_sol(gf_mat, BF, omega, C, K);
%% Residual
lhs     = del2(u_appr, dx, dy) + 1 ;
toc

%% Plot the numerical solution
figure(3),
surf(x, y, u_appr)
shading interp
%lighting phong
title('Numerical solution u_{appr}(x,y)','Interpreter','tex')
axis equal
%camlight('left');

% hold on
% stem3(x, y, u_appr)

figure(4),
surf(x(iC,iC), y(iC,iC), lhs(iC,iC))
shading interp
%lighting phong
% material metal 
title('Residual {res}(x,y)=\Delta u_{appr}(x,y)+1','Interpreter','tex')
axis square
%camlight('left');