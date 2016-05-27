/*			Petrov-Rvachev-Galerkin Method Solver			*/
/*			Basis: Cubic B-splines 					*/
/*			Author: V.V.Bondarenko					*/
/*			Date: 20.05.2016					*/
/*			Version: 1.0.2						*/


#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>

#define N 5
#define A -1.
#define B 1.
double right_part_f(double x, double y)
{
	return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));
/*	return 2.*(x*x+y*y-2.);*/
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

double integralRight(double (*f)(double, double, int), double x0, double x1, double y0, double y1, int m)
{
	double i,j, stepx=(x1-x0)/50., stepy=(y1-y0)/50.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{
			res += stepx*stepy*(*f)(i,j,m);
		}
	}
	return res;
}

double omega(double x, double y)
{
//	double result = (M_PI*M_PI-x*x)*(M_PI*M_PI-y*y);
	double result = (1.-x*x)*(1.-y*y);

/*	if(result < 0.) */
/*		return 0.;*/
	return result;
}

double phi(double x, double y, int n)
{
	return f_B_3(0.5*(x-(double)(n%N))/(B-A))*f_B_3(0.5*(y-(double)(n/N))/(B-A));
	//return pow(x,n%N)*pow(y,n/N);
}

double basis(double (*f)(double, double,int),double x, double y, int n)
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

void form_matrix(gsl_matrix * system, gsl_vector * RightPart, double x1, double x2, double y1, double y2)
{
	int i, j;
	for(i = 0; i < N*N; i++)
	{
		gsl_vector_set(RightPart, i, -integralRight(right_under_int,x1,x2,y1,y2,i));
		for(j = 0; j < N*N; j++)
		{
			gsl_matrix_set(system, i,j, integralLeft(left_under_int,x1,x2,y1,y2,i,j));
		}
	}
}

void solve_matrix_eq(gsl_vector * solution, gsl_matrix * system, gsl_vector * RightPart)
{
	int i;
	gsl_permutation * p = gsl_permutation_alloc (N*N);
	gsl_linalg_LU_decomp (system, p, &i);
	gsl_linalg_LU_solve (system, p, RightPart, solution);
}

double reconstruct_at(gsl_vector *solution, double x, double y)
{
	int i; double result = 0.;
	for(i=0; i<N*N; i++)
	{
		result += gsl_vector_get(solution, i)*basis(phi,x,y,i);
	}
	return result;
}

double plot_region(gsl_vector *solution, double x1, double x2, double y1, double y2)
{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("plot_region", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, reconstruct_at(solution,i,j));
	fclose(op);
}

int main()
{
	double a = A, b =B;
	gsl_matrix *sys = gsl_matrix_alloc (N*N,N*N);;
	gsl_vector  *rightpart	= gsl_vector_alloc(N*N),
			*solution	= gsl_vector_alloc(N*N);
	
	form_matrix(sys, rightpart, a,b, a,b);
	solve_matrix_eq(solution, sys, rightpart);
	
	plot_region(solution, a,b, a,b);
	
	return 0;
}
