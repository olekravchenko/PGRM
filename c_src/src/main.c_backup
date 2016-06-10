#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>
int N;
double intStep, glob_delta;

#include "right_parts.c"
#include "basis_functions.c"


double diff_step;



double Simpson_int(double (*f)(double,double, int, int), double x0, double x1, 
                                    double y_c, int k1, int k2)
{
	double 	i,
			step = (x1-x0)/intStep;
			//res = 0.;
      const int   q_of_layers = (int)((x1-x0)/step) + 1;
	int         k;
	double res = (*f)(x0, y_c, k1, k2) + 
	             (*f)((double)(q_of_layers-1)*step + x0, y_c, k1, k2);
	for(k = 1; k < q_of_layers-1; k++)
	{
            res += 2.*(double)(1+k%2)*(*f)((double)k*step + x0, y_c, k1, k2);
	}
	return step/3.*res;
}
//Simpson_int_int
double integralLeft(double (*f)(double,double, int, int),
                  double x0, double x1, double y0, double y1, int k1, int k2)
{
	double	i,
			step = (y1-y0)/intStep;
			//res = 0.;
	const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
	int         k;
	double res = Simpson_int((*f),x0,x1,y0, k1,k2) + 
	           Simpson_int((*f),x0,x1,y0+(double)(q_of_layers-1)*step, k1,k2);
	for(k = 1; k < q_of_layers-1; k++)
	{
            res += 2.*(double)(1+k%2)*Simpson_int((*f),x0,x1,y0+(double)k*step, k1,k2);
	}
	return step/3.*res;
}

double Simpson_right(double (*f)(double,double, int), double x0, double x1, 
                                    double y_c, int k1)
{
	double 	i,
			step = (x1-x0)/64.;
			//res = 0.;
      const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
	int         k;
	double res = (*f)(x0, y_c, k1) + 
	             (*f)((double)(q_of_layers-1)*step + x0, y_c, k1);
	for(k = 1; k < q_of_layers-1; k++)
	{
            res += 2.*(double)(1+k%2)*(*f)((double)k*step + x0, y_c, k1);
	}
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
			step = (y1-y0)/64.;
			//res = 0.;
	const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
	int         k;
	double res = Simpson_right((*f),x0,x1,y0, k1) + Simpson_right((*f),x0,x1,y0+(double)(q_of_layers-1)*step, k1);
	for(k = 1; k < q_of_layers-1; k++)
	{
            res += 2.*(double)(1+k%2)*Simpson_right((*f),x0,x1,y0+(double)k*step, k1);
	}
	return step/3.*res;
}







double omega(double x, double y)
// Returns value of R-function \omega(x,y)
// ToDo: modify for random bound positions
{
	return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
	//return result;
}

double basis(double (*f)(double, double,int),double x, double y, int n)
// Returns value of n-th \psi-basis function, used further, at point (x,y)
{
	return (*f)(x,y,n)*omega(x,y);
}


//testing variant
double left_under_int(double x, double y, int m, int n)
{

	double 	omega_px = omega(x + diff_step, y),
			omega_mx = omega(x - diff_step, y),
			omega_py = omega(x, y + diff_step),
			omega_my = omega(x, y - diff_step),
			phim_px  = phi(x+diff_step,y,m),
			phim_py  = phi(x,y+diff_step,m),
			phim_mx  = phi(x-diff_step,y,m),
			phim_my  = phi(x,y-diff_step,m),
			phin_px  = phi(x+diff_step,y,n),
			phin_py  = phi(x,y+diff_step,n),
			phin_mx  = phi(x-diff_step,y,n),
			phin_my  = phi(x,y-diff_step,n);
	
	return 0.25*glob_delta*glob_delta*
	(     (omega_px*phin_px - omega_mx*phin_mx)*
	      (omega_px*phim_px - omega_mx*phim_mx)+
	      (omega_py*phin_py - omega_my*phin_my)*
	      (omega_py*phim_py - omega_my*phim_my)
	);
}


double left_under_int_old(double x, double y, int m, int n)
// Returns \nabla\psi_m \nabla\psi_n for integral calculation
{
	double res, delta = diff_step;
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

void errors_to_stdio
     (gsl_vector *solution, 
      double x1, double x2, 
      double y1, double y2)
{
	double	hx = (x2-x1)/64.,
			hy = (y2-y1)/64.,
			i,j,maxerr=0.,err;

	for(i=x1; i<=x2; i+=hx)
		for(j=y1; j<=y2; j+=hy)
			if((err = fabs(reconstruct_at(solution,i,j)-u_exact(i,j)))>maxerr) maxerr = err;
	printf("%d %g %f\n",N*N,intStep,maxerr);
}

int main(int argc, char **argv)
{
	//double a = A, b = B;
	N=atoi(argv[1]);
	intStep = (double) atoi(argv[2]);
	
	init_eq(atoi(argv[3]));
	init_basis(atoi(argv[4]));
	diff_step = pow(2.,-9);
	glob_delta = 1./diff_step;
	
	gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
	gsl_vector  *rightpart	= gsl_vector_alloc(N*N),
			*solution	= gsl_vector_alloc(N*N);
	
	form_matrix		(sys, rightpart, X0,X1, Y0,Y1);
	
/*	FILE *op;*/
/*	op = fopen("./matrix", "w");*/
/*	gsl_matrix_fprintf(op, sys, "%f");*/
/*	fclose(op);*/
	
	solve_matrix_eq	(solution, sys, rightpart);
	
	errors_to_stdio	(solution, X0,X1, Y0,Y1);
//	plot_region		(solution, X0,X1, Y0,Y1);
//	plot_region_error	(solution, X0,X1, Y0,Y1);
//	plot_exact_solution	(X0,X1, Y0,Y1);
//	plot_omega			(X0,X1, Y0,Y1);

//	system("./Plot");
//	system("./Plot_err");
//	system("./Plot_exact");
//	system("./Plot_omega");
	return 0;
}
