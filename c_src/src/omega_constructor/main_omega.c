#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct omega_primitive{
	double x0, y0;
	double a, b;
	int type;
}omega_primitive;

/*
 * полоса
 * прямоугольник
 * круг (эллипс)
 * треугольник
 * полуплоскость
 * 
 * 
 * парабола
 * угол
 * 
 */
