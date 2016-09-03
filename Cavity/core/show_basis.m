function [] = show_basis(omega, nf, hf, x, y)

for i = 1:nf+2
    for j = 1:nf+2
        k = (i-1)*(nf+2) + j;
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