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
a = -pi; b = -a;
c = -pi; d = -c;

nx = 2^5;
ny = nx;
dx = (b-a)/(nx-1);
dy = (d-c)/(ny-1);

tx = a:dx:b;
ty = c:dy:d;
[x, y] = meshgrid(tx, ty);

%% Input functions
%rhs_func    = @(x,y) 2*x.*y.*exp(x+y).*(x+y+x.*y-3);    % Right hand side function
%ua          = @(x,y) x.*(x-1).*y.*(y-1).*exp(x+y);      % Exact Solution

rhs_func    = @(x,y) -2*sin(x).*sin(y);    % Right hand side function
ua          = @(x,y) sin(x).*sin(y);      % Exact Solution




%% Omega function
% om1 = 0.5*(0.5 - (x - 0.5).^2 - (y - 0.5).^2 - ...
%     sqrt(x.^2.*(x-1.0).^2 + y.^2.*(y-1.0).^2));
% om2 = 1.5 - x - y;
% om3 = y - 2*x + 1.5;
% om12 = 0.5*(om1 + om2 - sqrt(om1.^2 + om2.^2));
% omega = 0.5*(om12 + om3 - sqrt(om12.^2 + om3.^2));
% omega = (0.25 - (x - 0.5).^2) .* (0.25 - (y - 0.5).^2);
% omega(omega<0) = 0;
% omega = x.*(x-b).*y.*(y-d);
om1 = x + b; om2 = y + d; 
om3 = b - x; om4 = d - y;
om12 = r_con(om1, om2); 
om34 = r_con(om3, om4);
omega = r_con(om12, om34);

u_appr = 0*omega;

% om_ch = omega;
% om_ch(om_ch>0) = 1.0;

% figure(11),
% surf(x, y, om_ch)
% shading interp
% lighting phong
% title('Characteristic function of \omega(x,y)','Interpreter','latex')
% axis square
% % view(0,90)
% % axis off
% xlabel('x','Interpreter','latex')
% ylabel('y','Interpreter','latex')
% zlabel('\widetilde{\omega}(x,y)','Interpreter','latex')
% 
% u_appr = 0*omega;

figure(1),
surf(x, y, omega)
shading interp;  axis square; %axis off; 
title('\omega(x,y)','Interpreter','tex')
% axis equal

figure(2),
nf = 4;
for i = 1:nf
    for j = 1:nf
        k = (i-1)*nf + j;
        subplot(nf,nf,k)
        tmp = basis_function(k, x, y);
        surf(x, y, omega.*tmp)
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
K   = 16;
A   = zeros(K,K);
B   = zeros(1,K);
BF  = zeros(K,nx,ny);
ind = 1:K;

%% Computational Loop
for l = ind
   BF(l,:,:)   = basis_function(l, x, y);
   B(l)        = trapz(tx, trapz(ty, rhs_func(x,y).*omega.*squeeze(BF(l,:,:))));    
   for k = ind                
       tmp     = del2(omega.*squeeze(BF(l,:,:)), dx, dy);
       A(k, l) = trapz(tx, trapz(ty, omega.*squeeze(BF(k,:,:)).*tmp));
   end % for l
end % for k

% figure(2),
% spy(A)

C = cgs(A,B');

%% Reconstruction of numerical solution
for k = ind
    u_appr = u_appr + C(k)*squeeze(BF(k,:,:)).*omega.*omega./dx;
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
lhs = del2(u_appr, dx, dy) - rhs_func(x,y);

figure(5),
surf(x, y, lhs)
shading interp;  
title('Residual res(x,y)=\Delta u_{appr}(x,y)+1','Interpreter','tex')
axis square