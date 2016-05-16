function [] = show_basis(omega, nf, hf, x, y)

nx = length(x); ny = length(y);

for i = 1:nf+2
    for j = 1:nf+2
        k = (i-1)*(nf+2) + j;
        subplot(nf+2, nf+2, k)
        tmp01 = bdspln_function(i-1, (x-x(1,1))/(x(nx,nx)-x(1,1)), hf, nf);
        tmp02 = bdspln_function(nf-j+1, (y-y(1,1))/(y(ny,ny)-y(1,1)), hf, nf);
        tmp = tmp01 .* tmp02;
        surf(x, y, omega .* tmp)
        shading interp
%        lighting phong
%        camlight('left')
        axis square
        axis off
        view(0,90)
        xlabel('x','Interpreter','tex')
        ylabel('y','Interpreter','tex')
        title(['\psi_{' int2str(k) '}=' '\omega {B}_{3,' int2str(i-1) '}{B}_{3,' int2str(j-1) '}'],'Interpreter','tex')
%         title(['(' int2str(i-1),',' int2str(j-1) ')'])
    end
end
% text(-5,-0.5,'B-spline basis function $\mathrm{B}_{3,k}(x,y)=\mathrm{B}_{3,i}(x)\mathrm{B}_{3,j}(y)$','Interpreter','latex')