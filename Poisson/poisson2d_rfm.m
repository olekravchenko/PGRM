%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 2-D Poisson equation with
%           Petrov-Galerkin-Rvachev method (PGRM)
%
%               \nabla^2 u(x,y) = f(x,y),
%               u(x,y)|_{\partial\Omega}  = 0,              (Dirichlet problem)
%               where \Omega domain is defined by R(Rvachev)-functions
%           
%           Here u^N(x,y) = \sum_k^N \psi_k(x,y), \psi_k = x^iy^j, 
%               k is multiindex, i=i(k), j=j(k).
%
%              coded by Oleg Kravchenko, BMSTU, 2015.07.02
%              UPD1: 2016.05.16
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Petrov Galerkin Scheme Implementation with the usage of R(Rvachev)-functions method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Pathes
addpath('core')
addpath('core_rfm')

%% Initialization
a = -1; b = -a;
c = -1; d = -c;

nx = 2^6;
ny = nx;
dx = (b-a)/(nx-1);
dy = (d-c)/(ny-1);

tx = a:dx:b;
ty = c:dy:d;
[x, y] = meshgrid(tx, ty);

K   = 3;            % 
nf  = K^2;          % Number of basis functions

%% Input functions
id_func = 1; id_omega = 1;
[ua, rhs_func]  = GetFunc(id_func);
omega           = GetOmega(id_omega, b, d);
om              = omega(x,y);
u_appr          = 0*om;

figure(1),
surf(x, y, om)
shading interp;  axis square; %axis off; 
title('\omega(x,y)','Interpreter','tex')
% axis equal

figure(2),
% nf = 16;
for i = 1:K
    for j = 1:K
        k = (i-1)*K + j;
        subplot(K,K,k)
        tmp = basis_function(k, x, y);
        surf(x, y, om.*tmp)
        shading interp;  axis square; axis off; 
%         view(36,33)
        view(0,90)
        xlabel('x','Interpreter','tex')
        ylabel('y','Interpreter','tex')
        title(['\psi_{' int2str(k) '}=' '\omega x^{' int2str(i-1) '}y^{' int2str(j-1) '}'],'Interpreter','tex')
    end
end
%text(-5,-0.5,'Basis function \psi_k(x,y)=x^iy^j','Interpreter','latex')

%% Assembling
A   = zeros(K,K);
B   = zeros(1,K);
BF  = zeros(K,nx,ny);
ind = 1:K;

%% Computational Loop
for l = ind
   BF(l,:,:)   = basis_function(l, x, y);
   B(l)        = trapz(tx, trapz(ty, rhs_func(x,y).*om.*squeeze(BF(l,:,:))));    
   for k = ind                
       tmp     = 4*del2(om.*squeeze(BF(k,:,:)), dx, dy);
       A(k, l) = trapz(tx, trapz(ty, om.*squeeze(BF(l,:,:)).*tmp));
   end % for l
end % for k

% figure(2),
% spy(A)
% A = sparse(A);
% C = cgs(A', B');
C = A' \ B';

%% Reconstruction of numerical solution
for k = ind
    u_appr = u_appr + C(k)*squeeze(BF(k,:,:)).*om;
end % for k

%% Plot the numerical solution
figure(3),
surf(x, y, u_appr)
shading interp;  
title('Numerical solution u_{appr}(x,y)','Interpreter','tex')
axis square
xlabel('x','Interpreter','tex')
ylabel('y','Interpreter','tex')
zlabel('u(x,y)','Interpreter','tex')

figure(4)
surf(x, y, ua(x,y))
shading interp;  
title('Exact solution u(x,y)','Interpreter','tex')
axis square
xlabel('x','Interpreter','tex')
ylabel('y','Interpreter','tex')
zlabel('u(x,y)','Interpreter','tex')

%% Errors
lhs = abs(4*del2(u_appr, dx, dy) - rhs_func(x,y));

disp(['Absolute error: ' num2str(max(max(lhs)))])

figure(5),
surf(x, y, lhs)
shading interp;  
title('Residual res(x,y)=\Delta u_{appr}(x,y)-f(x,y)','Interpreter','tex')
axis square