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

#include "simpson_integrals.c"


double omega(double x, double y)
// Returns value of R-function \omega(x,y)
// ToDo: modify for random bound positions
{
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}

double basis(double x, double y, int n)
// Returns value of n-th \psi-basis function, used further, at point (x,y)
{
	return phi(x,y,n)*omega(x,y);
}
#include "plotters.c"

//testing variant
double left_under_int(double x, double y, int m, int n)
{

    double  omega_px = omega(x + diff_step, y),
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


void form_matrix(gsl_matrix * system,
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
        gsl_vector_set(RightPart, i, -integralRight(right_part_f,basis,x1,x2,y1,y2,i));
        for(j = 0; j < N*N; j++)
        {
            gsl_matrix_set(system, i,j, integralLeft(left_under_int,x1,x2,y1,y2,i,j));
        }
    }
}

void solve_matrix_eq(gsl_vector * solution,
                     gsl_matrix * system,
                     gsl_vector * RightPart)
//Solve SLE Ax=b, where A = system, b = RightPart, x = solution
{
    int i;

    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (system, p, &i);
    gsl_linalg_LU_solve (system, p, RightPart, solution);
}


void errors_to_stdio(gsl_vector *solution,
                     double x1, double x2,
                     double y1, double y2)
{
    double  hx = (x2-x1)/64.,
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

    return 0;
}
