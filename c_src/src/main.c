#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>

#define N 5

double right_part_f(double x, double y)
{
	return -2.*sin(x)*sin(y);
}

double traps(double (*f)(double), double x0, double x1) 
//Trapezoidal rule integration for 1D
{
	double i;
	double step =(x1-x0)/50.;
	double ret = (*f)(x0);
	for(i=x0+step;i<x1; i+=step)
	{
		ret += 2.*(*f)(i);
	}
	ret *=0.5*step;
	return ret;
}

double integralLeft(double (*f)(double, double, int, int), double x0, double x1, double y0, double y1, int k1, int k2)
{
	double i,j, stepx=(x1-x0)/50., stepy=(y1-y0)/50.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{
			res += stepx*stepy*(*f)(i,j,k1,k2);
		}
	}
	return res;
}

double integralRight(double (*f)(double, double), double x0, double x1, double y0, double y1)
{
	double i,j, stepx=(x1-x0)/50., stepy=(y1-y0)/50.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{
			res += stepx*stepy*(*f)(i,j);
		}
	}
	return res;
}

double omega(double x, double y)
{
	return (M_PI*M_PI-x*x)*(M_PI*M_PI-y*y);
}

double phi(double x, double y, int n)
{
	return pow(x,n%N)*pow(y,n/N);
}

double basis(double (*f)(double, double),double x, double y, int n)
{
	return (*f)(x,y,n)*omega(x,y);
}

double left_under_int(double x, double y, int m, int n)
{
	double res, delta = 0.01;
	//int i, j;
	
	res = 0.25/delta/delta*
(	(basis(phi, x+delta, y, n)-basis(phi, x-delta, y, n))*
	(basis(phi, x+delta, y, m)-basis(phi, x-delta, y, m))+
 	(basis(phi, x, y+delta, n)-basis(phi, x, y-delta, n))*
 	(basis(phi, x, y+delta, m)-basis(phi, x, y-delta, m))
);
	return res;
}

double right_under_int(double x, double y, int n)
{
	return right_part_f(x,y)*basis(phi, x,y, n);
}

void form_matrix(gsl_matrix * system, gsl_vector * RightPart double x1, double x2, double y1, double y2)
{
	int i, j;
	for(i = 0; i < N; i++)
	{
		gsl_vector_set(RightPart, integralRight());
		for(j = 0; j < N; j++)
		{
			
		}
	}

}

int main()
{
	
	return 0;
}
