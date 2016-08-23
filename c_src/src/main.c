#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_poly.c"
#include <string.h>
#include <pthread.h>
#include "smplDynArray.c"
#include <stdbool.h>
int N;
bool FiniteBasis = false;
double intStep, glob_delta, diff_step;

#include "tasks.c"
#include "basis_functions.c"
#include "plotters.c"
#include "gauss_integrals.c" //todo: reduce file to single function
#include "error_functions.c" 

double left_under_int(double x, double y, int m, int n)
{
 // \phi_m \Delta \phi_n 
    return  basis(x,y,m)*(
			basis(x+diff_step,y,n)+basis(x-diff_step,y,n)+
			basis(x,y+diff_step,n)+basis(x,y-diff_step,n)
			-4.*basis(x,y,n))*glob_delta*glob_delta;
}
double laplacian(double x, double y, int n)
{
 // \phi_m \Delta \phi_n 
    return (basis(x+diff_step,y,n)+basis(x-diff_step,y,n)+
			basis(x,y+diff_step,n)+basis(x,y-diff_step,n)
			-4.*basis(x,y,n))*glob_delta*glob_delta;
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
void form_matrix_colloc (gsl_matrix * system,
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
        gsl_vector_set(RightPart, i, right_part_f(X0+cubic_stepx*(double)(1+i/N),
												  Y0+cubic_stepy*(double)(1+i%N)));
        for(j = 0; j < N*N; j++)
        {
            gsl_matrix_set(system, i,j, laplacian(X0+cubic_stepx*(double)(1+i/N),
												  Y0+cubic_stepy*(double)(1+i%N),j));
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
    //form_matrix		(sys, rightpart, X0,X1, Y0,Y1);
    form_matrix_colloc		(sys, rightpart, X0,X1, Y0,Y1);
    //form_matrix_parallel(sys, rightpart, X0,X1, Y0,Y1);
	FILE* matr_op;
	matr_op = fopen("matrix.txt","w");
	gsl_matrix_fprintf(matr_op,sys,"%3.3f");
	fclose(matr_op);


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
