#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_poly.c"
#include <string.h>
#include <pthread.h>
#include "smplDynArray.c"
int N;
double intStep, glob_delta;

#include "right_parts.c"
#include "basis_functions.c"

double diff_step;


/*
double omega(double x, double y)
// Returns value of R-function \omega(x,y)
// ToDo: modify for random bound positions
{
	return (442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
     Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y) - 
     225*Sqrt(Power(-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y),2) + 
        Power(-8 + 225*x*x + 225*y*y + 
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y),2)/50625.))/225.;
}*/

double basis(double x, double y, int n)
// Returns value of n-th \psi-basis function, used further, at point (x,y)
{
		return phi(x,y,n)*omega(x,y);
}

#include "plotters.c"
#include "gauss_integrals.c" //todo: reduce file to single function
#include "error_functions.c" 


double left_under_int(double x, double y, int m, int n)
{
 // \phi_m \Delta \phi_n 
    double  omega_0  = omega(x,y),
			omega_px = omega(x + diff_step, y),
            omega_mx = omega(x - diff_step, y),
            omega_py = omega(x, y + diff_step),
            omega_my = omega(x, y - diff_step),
            phin_px  = phi(x+diff_step,y,n),
            phin_py  = phi(x,y+diff_step,n),
            phin_mx  = phi(x-diff_step,y,n),
            phin_my  = phi(x,y-diff_step,n);

    return omega_0*phi(x,y,m)*glob_delta*glob_delta*
			(omega_px*phin_px + omega_mx*phin_mx + 
			omega_py*phin_py + omega_my*phin_my - 4.*omega_0*phi(x,y,n));
    
}

/*
 * 
 * Place for Including of parallel former, if required
 * 
 */

void form_matrix (gsl_matrix * system,
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
        gsl_vector_set(RightPart, i, integralRight(right_part_f,basis,x1,x2,y1,y2,i));
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



int main(int argc, char **argv)
/*
 * requires 4 arguments to launch:
 * 
 * 1. N - quantity of basis functions per side. 
 * 2. intStep - quantity of integration nodes.
 * 3. id of example equation in right_parts.c.
 * 4. id of basis functions.
 * 
 * Typical launch command (one of the best for testing):
 * 	./bin/main 8 64 3 3
 * 
 * To build a program completly:
 * [path_to_PGRM/c_src]/build
 * or just double click on build script the same way as usual program.
 */
{
    //double a = A, b = B;
    N = atoi(argv[1]);
    intStep = (double) atoi(argv[2]);
	initGaussInt();
    init_eq(atoi(argv[3]));
    init_basis(atoi(argv[4]));
    diff_step = pow(2.,-9);
    glob_delta = 1./diff_step;

    gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
    gsl_vector  *rightpart	= gsl_vector_alloc(N*N);
    //*solution	= gsl_vector_alloc(N*N);
    solution_glob = gsl_vector_alloc(N*N);
    form_matrix		(sys, rightpart, X0,X1, Y0,Y1);
    //form_matrix_parallel(sys, rightpart, X0,X1, Y0,Y1);

    solve_matrix_eq	(solution_glob, sys, rightpart);
    //solution_glob = solution;
    //errors_to_stdio	(solution_glob, X0,X1, Y0,Y1);
    //multiplot		(solution_glob, X0,X1, Y0,Y1);
    plot_region		(solution_glob, X0,X1, Y0,Y1);
    //plot_region_error	(solution_glob, X0,X1, Y0,Y1);
    //plot_exact_solution	(X0,X1, Y0,Y1);
    //plot_omega		(X0,X1, Y0,Y1);

    return 0;
}
