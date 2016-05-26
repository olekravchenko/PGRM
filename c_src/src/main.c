#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>

#define N 6
#include "right_parts.c"
#include "basis_functions.c"






double Simpson_int(double (*f)(double,double, int, int), double x0, double x1, 
                                    double y_c, int k1, int k2)
{
	double 	i,
			step = (x1-x0)/512.,
			res = 0.;
      const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
	int         k;
	double      integrals_by_layer[q_of_layers];
	
	for(k = 0; k < q_of_layers; k++)
	{
	      integrals_by_layer[k] = (*f)((double)k*step + x0, y_c, k1, k2);
	}
	
	res = integrals_by_layer[0];
	for(k = 1; k < q_of_layers-1; k++)
	{
            if(k%2 == 0)
                  res += 2.*integrals_by_layer[k];
            else
                  res += 4.*integrals_by_layer[k];
	}
	res +=integrals_by_layer[q_of_layers-1];
	return step/3.*res;
}
//Simpson_int_int
double integralLeft(double (*f)(double,double, int, int),
                  double x0, double x1, double y0, double y1, int k1, int k2)
{
	double	i,
			step = (y1-y0)/512.,
			res = 0.;
	const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
	int         k;
	double      integrals_by_layer[q_of_layers];
	
	for(k = 0; k < q_of_layers; k++)
	{
	      integrals_by_layer[k] = Simpson_int((*f),x0,x1,y0+(double)k*step, k1,k2);
	}
	
	res = integrals_by_layer[0];
	for(k = 1; k < q_of_layers-1; k++)
	{
            if(k%2 == 0)
                  res += 2.*integrals_by_layer[k];
            else
                  res += 4.*integrals_by_layer[k];
	}
	res +=integrals_by_layer[q_of_layers-1];
	return step/3.*res;
}

double Simpson_right(double (*f)(double,double, int), double x0, double x1, 
                                    double y_c, int k1)
{
	double 	i,
			step = (x1-x0)/128.,
			res = 0.;
      const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
	int         k;
	double      integrals_by_layer[q_of_layers];
	
	for(k = 0; k < q_of_layers; k++)
	{
	      integrals_by_layer[k] = (*f)((double)k*step + x0, y_c, k1);
	}
	
	res = integrals_by_layer[0];
	for(k = 1; k < q_of_layers-1; k++)
	{
            if(k%2 == 0)
                  res += 2.*integrals_by_layer[k];
            else
                  res += 4.*integrals_by_layer[k];
	}
	res +=integrals_by_layer[q_of_layers-1];
	return step/3.*res;
}


//simpson right integral
double integralRight
     (double (*f)(double, double, int), 
      double x0, double x1, 
      double y0, double y1, 
      int k1)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
	double	i,
			step = (y1-y0)/128.,
			res = 0.;
	const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
	int         k;
	double      integrals_by_layer[q_of_layers];
	
	for(k = 0; k < q_of_layers; k++)
	{
	      integrals_by_layer[k] = Simpson_right((*f),x0,x1,y0+(double)k*step, k1);
	}
	
	res = integrals_by_layer[0];
	for(k = 1; k < q_of_layers-1; k++)
	{
            if(k%2 == 0)
                  res += 2.*integrals_by_layer[k];
            else
                  res += 4.*integrals_by_layer[k];
	}
	res +=integrals_by_layer[q_of_layers-1];
	return step/3.*res;
}







double omega(double x, double y)
// Returns value of R-function \omega(x,y)
// ToDo: modify for random bound positions
{
	double result = (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
      
	//if(result <= 0.) 
	//	return 0.;
	return result;
}

double basis(double (*f)(double, double,int),double x, double y, int n)
// Returns value of n-th \psi-basis function, used further, at point (x,y)
{
	return (*f)(x,y,n)*omega(x,y);
}

double left_under_int(double x, double y, int m, int n)
// Returns \nabla\psi_m \nabla\psi_n for integral calculation
{
	double res, delta = 0.0000000001;
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
// Return f(x,y)\psi_n (x,y)
{
	return right_part_f(x,y)*basis(phi, x,y, n);
}

void form_matrix
     (gsl_matrix * system, 
      gsl_vector * RightPart, 
      double x1, double x2, 
      double y1, double y2)
// Forms SLE system
// system 	- left part matrix form of system
// RightPart- right part vector of coefficients
// x1, x2	- sizes of rectangle by x
// y1, y2	- sizes of rectangle by y
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

void solve_matrix_eq
     (gsl_vector * solution, 
      gsl_matrix * system, 
      gsl_vector * RightPart)
//Solve SLE Ax=b, where A = system, b = RightPart, x = solution
{
	int i;

	gsl_permutation * p = gsl_permutation_alloc (N*N);
	gsl_linalg_LU_decomp (system, p, &i);
	gsl_linalg_LU_solve (system, p, RightPart, solution);
}

double reconstruct_at(gsl_vector *solution, double x, double y)
// Reconstucts value of solution at point (x,y)
{
	int i; double result = 0.;
	for(i=0; i<N*N; i++)
	{
		result += gsl_vector_get(solution, i)*basis(phi,x,y,i);
	}
	return result;
}

double plot_region
     (gsl_vector *solution, 
      double x1, double x2, 
      double y1, double y2)
// Plot solution in rectangle region 
// from x1 till x2 by x, and from y1 till y2 by y

{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_region", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, reconstruct_at(solution,i,j));
	fclose(op);
}

double plot_exact_solution
     (double x1, double x2, 
      double y1, double y2)
// Plot solution in rectangle region 
// from x1 till x2 by x, and from y1 till y2 by y

{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_exact_solution", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, u_exact(i,j));
	fclose(op);
}

double plot_omega
     (double x1, double x2, 
      double y1, double y2)
// Plot solution in rectangle region 
// from x1 till x2 by x, and from y1 till y2 by y

{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_plot_omega", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, omega(i,j));
	fclose(op);
}

double plot_region_error
     (gsl_vector *solution, 
      double x1, double x2, 
      double y1, double y2)
// Plot abs error of solution in rectangle region 
// from x1 till x2 by x, and from y1 till y2 by y
{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j;
	FILE * op;
	op = fopen("../plot_data/plot_region_error", "w");
	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			fprintf(op, "%f %f %f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
	fclose(op);
}

int main()
{
	//double a = A, b = B;
	init_eq(1);
	init_basis(1);
	
	gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
	gsl_vector  *rightpart	= gsl_vector_alloc(N*N),
			*solution	= gsl_vector_alloc(N*N);
	
	form_matrix		(sys, rightpart, X0,X1, Y0,Y1);
	
	FILE *op;
	op = fopen("./matrix", "w");
	gsl_matrix_fprintf(op, sys, "%f");
	fclose(op);
	
	
	solve_matrix_eq	(solution, sys, rightpart);
	
	plot_region		(solution, X0,X1, Y0,Y1);
	plot_region_error	(solution, X0,X1, Y0,Y1);
	plot_exact_solution	(X0,X1, Y0,Y1);
	plot_omega			(X0,X1, Y0,Y1);

	system("./Plot");
	system("./Plot_err");
//	system("./Plot_exact");
//	system("./Plot_omega");
	return 0;
}
