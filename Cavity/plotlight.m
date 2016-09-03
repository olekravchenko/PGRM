function plotlight(x, y, u)

figure(1),
surf(x, y, u); view(30, 30);
% daspect([1,1,1])
shading interp
light
lighting phong
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Surface','Interpreter','latex')
