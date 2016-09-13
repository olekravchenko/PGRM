#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct omega_primitive {
    double x0, y0;
    double a, b;
} omega_primitive;

/*
 * полоса || Ox		band_y			+
 * полоса || Oy		band_x			+
 * прямоугольник	rectangle		+
 * круг 			circle
 * эллипс			ellipse			+
 * треугольник		triangle		-
 * полуплоскость	half_plane		+
 *
 *
 * парабола			parabola		+
 * угол				angle			-
 *
 */

double band_y(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double width = properties.a;

    return 0.5/width*(width-X)*(width+X);
}

double band_x(omega_primitive properties, double x, double y) {
    double Y = y - properties.y0;
    double height = properties.b;

    return 0.5/height*(height-Y)*(height+Y);
}

double dx_band_x(omega_primitive properties, double x, double y) {

    return 0;
}

double dy_band_x(omega_primitive properties, double x, double y) {
    double Y = y - properties.y0;
    double height = properties.b;

    return -1/height*Y;
}

double d2x_band_x(omega_primitive properties, double x, double y) {
    
    return 0;
}

double dxdy_band_x(omega_primitive properties, double x, double y) {
    
    return 0;
}

double d2y_band_x(omega_primitive properties, double x, double y) {
    double Y = y - properties.y0;
    double height = properties.b;

    return -1/height;
}

double rectangle(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double width = properties.a;
    double height = properties.b;

    double val01 = 0.5/width*(width-X)*(width+X);
    double val02 = 0.5/height*(height-Y)*(height+Y);


    return val01 + val02 - sqrt(val01*val01 + val02*val02);
}

/*
double dx_rectangle(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    
    double val01 = X;
    double val02 = sqrt(X*X+Y*Y);


    return 1 - val01/val02;
}
*/

double circle(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double radius = properties.a;

    return 0.5/radius*(radius*radius-X*X-Y*Y);
}

double dx_circle(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double radius = properties.a;

    return -X/radius;
}

double dy_circle(omega_primitive properties, double x, double y) {
    double Y = y - properties.y0;
    double radius = properties.a;

    return -Y/radius;
}

double d2x_circle(omega_primitive properties, double x, double y) {
    double radius = properties.a;

    return -1/radius;
}

double dxy_circle(omega_primitive properties, double x, double y) {
    double radius = properties.a;

    return 0;
}

double d2y_circle(omega_primitive properties, double x, double y) {
    double radius = properties.a;

    return -1/radius;
}

double ellipse(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double b2 = b*b;
    double X2 = X*X;
    double Y2 = Y*Y;

    double val01 = a2*b2 - b2*X2 - a2*Y2;
    double val02 = 0.5/(a*b)*sqrt(a2 + b2 - X2 - Y2);

    return val01/val02;
}


double half_plane(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double b2 = b*b;

    return (a*X+b*Y)/sqrt(a2+b2);
}

double dx_half_plane(omega_primitive properties, double x, double y) {
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double b2 = b*b;

    return a/sqrt(a2+b2);
}

double dy_half_plane(omega_primitive properties, double x, double y) {
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double b2 = b*b;

    return b/sqrt(a2+b2);
}

double d2x_half_plane(omega_primitive properties, double x, double y) {
   
    return 0;
}

double dxy_half_plane(omega_primitive properties, double x, double y) {   

    return 0;
}

double d2y_half_plane(omega_primitive properties, double x, double y) {
   
    return 0;
}

double parabola(omega_primitive properties, double x, double y) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double X2 = X*X;

    return (b-a*X2-Y)/sqrt(1+a2*X2);
}
