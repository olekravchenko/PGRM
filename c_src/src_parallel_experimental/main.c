#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>
#include <pthread.h>
#include "smplDynArray.c"
int N;
double intStep, glob_delta;

#include "right_parts.c"
#include "basis_functions.c"

double diff_step;



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
#include "simpson_integrals.c"
#include "error_functions.c"
//testing variant
double left_under_int_old(double x, double y, int m, int n)
{//classical weak formulation
 // \Nabla \phi_m \Nabla \phi_n
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
                 (omega_py*phim_py - omega_my*phim_my));
}

double left_under_int(double x, double y, int m, int n)
{
 // \phi_m \Delta \phi_n 
    double  omega_0  = omega(x,y),
			omega_px = omega(x + diff_step, y),
            omega_mx = omega(x - diff_step, y),
            omega_py = omega(x, y + diff_step),
            omega_my = omega(x, y - diff_step),
            /*phim_px  = phi(x+diff_step,y,m),
            phim_py  = phi(x,y+diff_step,m),
            phim_mx  = phi(x-diff_step,y,m),
            phim_my  = phi(x,y-diff_step,m);*/
            phin_px  = phi(x+diff_step,y,n),
            phin_py  = phi(x,y+diff_step,n),
            phin_mx  = phi(x-diff_step,y,n),
            phin_my  = phi(x,y-diff_step,n);

    /*return 0.25*glob_delta*glob_delta*
           (     (omega_px*phin_px - omega_mx*phin_mx)*
                 (omega_px*phim_px - omega_mx*phim_mx)+
                 (omega_py*phin_py - omega_my*phin_my)*
                 (omega_py*phim_py - omega_my*phim_my));
    */
    return omega_0*phi(x,y,m)*glob_delta*glob_delta*
			(omega_px*phin_px + omega_mx*phin_mx + 
			omega_py*phin_py + omega_my*phin_my - 4.*omega_0*phi(x,y,n));
    
}



typedef struct rectangle_area {
    double x1,x2;
    double y1,y2;
} rectangle_area;

typedef struct parallel_arg {
    rectangle_area area;
    int i, j;
    double result;
} parallel_arg;

typedef struct parallel_arg2 {
    rectangle_area area;
    int threadId;
    double **result;
} parallel_arg2;

void *parallel_pre_former(void *arg)
{
    parallel_arg *a = (parallel_arg *)arg;
    a->result = integralLeft(left_under_int,a->area.x1,a->area.x2,a->area.y1,a->area.y2,a->i,a->j);
    return 0;
}

void *parallel_pre_former2(void *arg)
{
    parallel_arg2 *a = (parallel_arg2 *)arg;
    int id = a->threadId;
    int i, j;
    for(i = 0; i<N*N/16; i++)
    {
        for(j = 0; j<N*N; j++)
            a->result[i+id*N*N/16][j] = integralLeft(left_under_int,a->area.x1,a->area.x2,a->area.y1,a->area.y2,i+id*N*N/16,j);
    }
    //a->result = integralLeft(left_under_int,a->area.x1,a->area.x2,a->area.y1,a->area.y2,a->i,a->j);
    return 0;
}

void form_matrix_parallel(gsl_matrix * system,
                          gsl_vector * RightPart,
                          double x1, double x2,
                          double y1, double y2)
// Forms SLE system
// system 	- left part matrix form of system
// RightPart- right part vector of coefficients
// x1, x2	- sizes of rectangle by x
// y1, y2	- sizes of rectangle by y
{
    int i, j, k, return_code;
    parallel_arg2 arg[16];
    rectangle_area general_area;
    general_area.x1 = x1;
    general_area.x2 = x2;
    general_area.y1 = y1;
    general_area.y2 = y2;

    double **prematrixarray;
    init2DArr(&prematrixarray, N*N, N*N);
    pthread_t threads[16];

    for(i = 0; i < 16; i++)
    {
        arg[i].area = general_area;
        arg[i].threadId = i;
        arg[i].result = prematrixarray;
        return_code=pthread_create(&threads[i],NULL,parallel_pre_former2, (void*)&arg[i]);
    }

    for(k = 0; k < 16; k++)
    {
        return_code = pthread_join(threads[k], NULL);

        for(i = 0; i<N*N/16; i++)
            for(j = 0; j < N*N; j++)
                prematrixarray[i+k*N*N/16][j]=arg[i].result[i+k*N*N/16][j];
    }
    return_code++;
    for(i = 0; i < N*N; i++)
    {
		//for non variational solver remove minus
        gsl_vector_set(RightPart, i, integralRight(right_part_f,basis,x1,x2,y1,y2,i));
        for(j = 0; j < N*N; j++)
            gsl_matrix_set(system, i,j, prematrixarray[i][j]);
    }
    free2DArr(&prematrixarray, N*N);
}
void form_matrix_parallel_v1_with_bug(gsl_matrix * system,
                                      gsl_vector * RightPart,
                                      double x1, double x2,
                                      double y1, double y2)
// Forms SLE system
// system 	- left part matrix form of system
// RightPart- right part vector of coefficients
// x1, x2	- sizes of rectangle by x
// y1, y2	- sizes of rectangle by y
{
    int i, j, return_code;
    parallel_arg arg[N*N][N*N];
    rectangle_area general_area;
    general_area.x1 = x1;
    general_area.x2 = x2;
    general_area.y1 = y1;
    general_area.y2 = y2;

    double **prematrixarray;
    init2DArr(&prematrixarray, N*N, N*N);
    pthread_t threads[N*N][N*N];
    for(i = 0; i < N*N; i++)
        for(j = 0; j < N*N; j++)
        {
            arg[i][j].area = general_area;
            arg[i][j].i = i;
            arg[i][j].j = j;
            return_code = pthread_create(&threads[i][j],NULL,parallel_pre_former, (void*)&arg[i][j]);
        }

    for(i = 0; i < N*N; i++)
        for(j = 0; j < N*N; j++)
        {
            return_code = pthread_join(threads[i][j], NULL);
            prematrixarray[i][j] = arg[i][j].result;
        }
    return_code++;
    for(i = 0; i < N*N; i++)
    {
        gsl_vector_set(RightPart, i, -integralRight(right_part_f,basis,x1,x2,y1,y2,i));
        for(j = 0; j < N*N; j++)
        {
            //gsl_matrix_set(system, i,j, integralLeft(left_under_int,x1,x2,y1,y2,i,j));
            gsl_matrix_set(system, i,j, prematrixarray[i][j]);
        }
    }
    free2DArr(&prematrixarray, N*N);
}

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

    init_eq(atoi(argv[3]));
    init_basis(atoi(argv[4]));
    diff_step = pow(2.,-9);
    glob_delta = 1./diff_step;

    gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
    gsl_vector  *rightpart	= gsl_vector_alloc(N*N);
    //*solution	= gsl_vector_alloc(N*N);
    solution_glob = gsl_vector_alloc(N*N);
    //form_matrix		(sys, rightpart, X0,X1, Y0,Y1);
    form_matrix_parallel(sys, rightpart, X0,X1, Y0,Y1);

    //FILE *op;
    //op = fopen("./matrix", "w");
    //gsl_matrix_fprintf(op, sys, "%f");
    //fclose(op);

    solve_matrix_eq	(solution_glob, sys, rightpart);
    //solution_glob = solution;
    errors_to_stdio	(solution_glob, X0,X1, Y0,Y1);
    //plot_region		(solution_glob, X0,X1, Y0,Y1);
    //plot_region_error	(solution_glob, X0,X1, Y0,Y1);
    //plot_exact_solution	(X0,X1, Y0,Y1);
    //plot_omega		(X0,X1, Y0,Y1);

    return 0;
}
