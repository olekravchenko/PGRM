%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Solving 2-D Poisson equation with
%           Petrov-Galerkin-Rvachev method (PGRM)
%
%               \nabla^2 u(x,y) = f(x,y),
%               u(x,y)|_{\partial\Omega}  = 0,              (Dirichlet problem)
%               where \Omega domain is defined by R(Rvachev)-functions
%           
%           Here u^N(x,y) = \sum_k^N \psi_k(x,y), \psi_k = B_i(x)B_j(y), 
%               k is multiindex, i=i(k), j=j(k),
%               B_i(x)B_j(y) is a tensor product of 1D cubic B-splines,
%               from [4] .
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
% [4] S. Beck. Cubic spline finite element method for solving Poisson's equation
%     on a square. M.Sc. Thesis. Colorado School of Mines. 2008/
%     link: https://dspace.library.colostate.edu/bitstream/handle/11124/353/Beck_mines_0052N_10391.pdf?sequence=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Petrov Galerkin Scheme Implementation with the usage of R(Rvachev)-functions method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

%% Pathes
addpath('core')
addpath('core_rfm')

%% Initialization
a = 0; b = 1;
c = 0; d = 1;

nx = 2^7;
ny = nx;
dx = (b-a) / (nx-1);
dy = (d-c) / (ny-1);

tx = a:dx:b;
ty = c:dy:d;
[x, y] = meshgrid(tx, ty);


%% Omega function
om1 = 0.5*(0.5 - (x - 0.5).^2 - (y - 0.5).^2 - ...
    sqrt(x.^2.*(x-1.0).^2 + y.^2.*(y-1.0).^2));
om2 = 1.5 - x - y;
om3 = y - 2*x + 1.5;
om12 = 0.5*(om1 + om2 - sqrt(om1.^2 + om2.^2));
omega = 0.5*(om12 + om3 - sqrt(om12.^2 + om3.^2));
% omega = (0.25 - (x - 0.5).^2) .* (0.25 - (y - 0.5).^2);
omega(omega<0) = 0;

omega = tanh(300*omega);

u_appr = 0*omega;

figure(1),
surf(x, y, omega)
shading interp
lighting phong
camlight('left');
title('$\omega_w(x,y)=\tanh(300\omega(x,y))$','Interpreter','latex')
axis equal

figure(2),
nf = 2;
hf = 1 / (nf + 1);
ind_mat = zeros(3,(nf+2)^2);

tic
for i = 1:nf+2
    for j = 1:nf+2
        k = (i-1)*(nf+2) + j;
        % fill matrix of indexes
        ind_mat(1,k) = k;
        ind_mat(2,k) = i;
        ind_mat(3,k) = j;
        subplot(nf+2, nf+2, k)
        tmp = bdspln_function(i-1, x, hf, nf).*bdspln_function(j-1, y, hf, nf);
        surf(x, y, omega .* tmp)
        shading interp
        lighting phong
        axis square
        axis off
        view(0,90)
        xlabel('$x$','Interpreter','latex')
        ylabel('$y$','Interpreter','latex')
        title(['$\psi_{' int2str(k) '}=' '\omega \mathrm{B}_{3,' int2str(i-1) '}\mathrm{B}_{3,' int2str(j-1) '}$'],'Interpreter','latex')
%         title(['(' int2str(i-1),',' int2str(j-1) ')'])
    end
end
% text(-5,-0.5,'B-spline basis function $\mathrm{B}_{3,k}(x,y)=\mathrm{B}_{3,i}(x)\mathrm{B}_{3,j}(y)$','Interpreter','latex')

%% Assembling
K   = (nf+2)^2;
A   = zeros(K,K);
B   = zeros(1,K);
% C   = B;
BF  = zeros(K,nx,ny);
ind = 1:K;

%% Compute matrix of basis functions
for k = ind
    i = ind_mat(2,k);
    j = ind_mat(3,k);
    BF(k,:,:)   = bdspln_function(i-1, x, hf, nf).*bdspln_function(j-1, y, hf, nf);
end

for l = ind
    B(l)        = - trapz(tx, trapz(ty, omega.*squeeze(BF(l,:,:))));    
    for k = ind                
        tmp     = del2(omega.*squeeze(BF(l,:,:)), dx, dy);
        A(k, l) =   trapz(tx, trapz(ty, omega.*squeeze(BF(k,:,:)).*tmp));
    end % for l
end % for k

% figure(2),
% spy(A)

C = cgs(A,B');

%% Reconstruction of numerical solution
for k = ind
    u_appr = u_appr + C(k)*squeeze(BF(k,:,:)).*omega;
end % for k

tic
%% Plot the numerical solution
figure(3),
surf(x, y, u_appr)
shading interp
lighting phong
title('Numerical solution $u_{appr}(x,y)$','Interpreter','latex')
axis square

%% Errors
lhs = del2(u_appr, dx, dy) + 1 ;

iC = 1:nx;

figure(4),
surf(x(iC,iC), y(iC,iC), lhs(iC,iC))
shading interp
lighting phong
title('Residual $res(x,y)=\triangle u_{appr}(x,y)+1$','Interpreter','latex')
axis square