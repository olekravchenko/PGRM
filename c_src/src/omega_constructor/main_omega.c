#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct omega_primitive {
    double x0, y0;
    double a, b;
} omega_primitive;


double R_and(double x, double y)
{
	return x+y-sqrt(x*x+y*y);
}
double R_or(double x, double y)
{
	return x+y+sqrt(x*x+y*y);
}

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


double half_plane(omega_primitive properties, double x, double y, int dxN, int dyN) 
{
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double a = properties.a;
    double b = properties.b;
	
	if(dxN == 0 && dyN == 0)
		return (a*X+b*Y)/sqrt(a*a+b*b);
	if(dxN == 1 && dyN == 0)
		return a/sqrt(a*a+b*b);
	if(dxN == 0 && dyN == 1)
		return b/sqrt(a*a+b*b);
		
	return 0.;
}

double band_y(omega_primitive properties, double x, double y, int dxN, int dyN) {
    double X = x - properties.x0;
    double width = properties.a;

	if(dxN == 0 && dyN == 0)
		return 0.5/width*(width-X)*(width+X);
	if(dyN > 0 || dxN > 2)
		return 0.;
	if(dxN == 1)
		return -X/width;
	if(dxN == 2)
		return -1./width;
	return 0.;
}

double band_x(omega_primitive properties, double x, double y, int dxN, int dyN) {
    double Y = y - properties.y0;
    double height = properties.b;
	
	if(dxN == 0 && dyN == 0)
		return 0.5/height*(height-Y)*(height+Y);
	if(dyN > 2 || dxN > 0)
		return 0.;
	if(dyN == 1)
		return -Y/height;
	if(dyN == 2)
		return -1./height;
	return 0.;
}



double circle(omega_primitive properties, double x, double y, int dxN, int dyN) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double radius = properties.a;


	if(dxN == 0 && dyN == 0)
		return 0.5/radius*(radius*radius-X*X-Y*Y);
	if(dxN == 0)
	{
		if(dyN == 1)
			return -Y/radius;
		if(dyN == 2)
			return -1./radius;
	}
	if(dyN == 0)
	{
		if(dxN == 1)
			return -X/radius;
		if(dxN == 2)
			return -1./radius;
	}
	return 0.;
}


double ellipse(omega_primitive properties, double x, double y, int dxN, int dyN) {
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







double rectangle(omega_primitive properties, double x, double y, int dxN, int dyN) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double width = properties.a;
    double height = properties.b;

    double val01 = 0.5/width*(width-X)*(width+X);
    double val02 = 0.5/height*(height-Y)*(height+Y);


    return val01 + val02 - sqrt(val01*val01 + val02*val02);
}

double line_segment(omega_primitive properties, double x, double y, int dxN, int dyN)
{
	double	A = 1./(properties.a-properties.x0),
			B = 1./(properties.b-properties.y0),
			norm = 1./sqrt(A*A+B*B);
	
	return norm*R_or(fabs(A*(x-properties.x0)-B*(y-properties.y0)),
			R_or(-A*(x-properties.x0)-B*(y-properties.y0),
					A*(x-properties.a)+B*(y-properties.b)));
	
}

double parabola(omega_primitive properties, double x, double y, int dxN, int dyN) {
    double X = x - properties.x0;
    double Y = y - properties.y0;
    double a = properties.a;
    double b = properties.b;

    double a2 = a*a;
    double X2 = X*X;

    return (b-a*X2-Y)/sqrt(1+a2*X2);
}



