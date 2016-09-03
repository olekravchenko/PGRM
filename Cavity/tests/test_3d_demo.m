close all
clear all
clc

[x,y,z,v] = flow;
isovalue = -1;
purple = [1.0 0.5 1.0];
figure;
p = patch(isosurface(x,y,z,v,isovalue));
isonormals(x,y,z,v,p);
set(p,'FaceColor',purple,'EdgeColor','none');
view([-10 40]);
axis on;
grid on;
light;
lighting phong;
camlight('left');