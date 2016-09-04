#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include <openacc.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_poly.c"
#include <string.h>
#include <pthread.h>
#include "smplDynArray.c"

typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;

typedef struct basis_args {
	double	x, y;
	int		m, n;
} basis_args;

//Definitions moved to tasks.h
//ToDo: reduce usage of global variables
#include "tasks.h"
#include "basis_functions.c"
#include "plotters_new.c"
#include "gauss_integrals_compact.c"


//#include "error_functions.c" 
//ToDo: avoid usage of globally defined
//		re-write code in error_functions.c

double left_under_int_new(basis_args arguments)
{
	double 	x = arguments.x;
	double 	y = arguments.y;
	int 	m = arguments.m;
	int 	n = arguments.n;
	
    return  structure(x,y,m)*(
			structure(x+diff_step,y,n)+structure(x-diff_step,y,n)+
			structure(x,y+diff_step,n)+structure(x,y-diff_step,n)
			-4.*structure(x,y,n))*glob_delta*glob_delta;
}

double right_under_int_new(basis_args arguments)
{
	double 	x = arguments.x;
	double 	y = arguments.y;
	int 	m = arguments.m;
	
	return right_part_f(x,y)*structure(x,y,m);
}

void form_matrix_new (gsl_matrix * system,
                      gsl_vector * RightPart,
                      rect_area int_area)
// Forms SLE system
// system 	- left part matrix form of system
// RightPart- right part vector of coefficients
// x1, x2	- sizes of rectangle by x
// y1, y2	- sizes of rectangle by y
{
    int i, j;

    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
#pragma omp parallel for shared(int_area, system, RightPart,N) private(args, i,j)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(RightPart, i, gauss_integral(right_under_int_new,int_area,args,2));
        for(j = 0; j < N*N; j++)
        {
            args.n = j;
            gsl_matrix_set(system, i,j, gauss_integral(left_under_int_new, int_area, args,2));
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
 * Typical launch command with arguments:
 * 	./main 8 4 3 3
 * 
 * --//-- w/o arguments:
 *  ./main
 * 
 * Launch w/o arguments equal to:
 *  ./main 8 4 10 5
 * 
 * To build a program completly:
 * [path_to_PGRM/c_src]/src/build.sh
 * or just double click on build script the same way as usual program.
 * 
 * ToDo: 
 * 	-add normal argument processing
 */
{
    initGaussInt();

    omp_set_dynamic(1);      
    omp_set_num_threads(16);


    if(argc>=5)
    {
        N 			= atoi(argv[1]);
        intStep 	= (double) atoi(argv[2]);
        init_eq(atoi(argv[3]));
        init_basis(atoi(argv[4]));
    }
    else if(argc == 1)
    {
        N 			= 8;
        intStep 	= 4.;
        init_eq(10);
        init_basis(5);
	}
    
    diff_step 	= pow(2.,-9);
    glob_delta 	= 1./diff_step;
    
	rect_area sol_area;
	sol_area.x0 = X0;
	sol_area.x1 = X1;
	sol_area.y0 = Y0;
	sol_area.y1 = Y1;
	
	
    gsl_matrix 	*sys 		= gsl_matrix_alloc (N*N,N*N);;
    gsl_vector  *rightpart	= gsl_vector_alloc(N*N),
				*solution	= gsl_vector_alloc(N*N);
	
	
    form_matrix_new	(sys, rightpart, sol_area);
    //form_matrix			(sys, rightpart, X0,X1, Y0,Y1);
    
	//FILE* matr_op;
	//matr_op = fopen("matrix.txt","w");
	//gsl_matrix_fprintf(matr_op,sys,"%3.3f");
	//fclose(matr_op);


    solve_matrix_eq	(solution, sys, rightpart);
    //errors_to_stdio	(solution_glob, X0,X1, Y0,Y1);
    //multiplot		(solution_glob, X0,X1, Y0,Y1);
    
    plot_region_colorplot(solution, sol_area);
    //plot_region_colorplot(solution, X0,X1, Y0,Y1);
    //plot_laplacian	(solution, X0,X1, Y0,Y1);
    //plot_region_error	(solution_glob, X0,X1, Y0,Y1);
    //plot_exact_solution	(X0,X1, Y0,Y1);
    //plot_omega		(X0,X1, Y0,Y1);

    return 0;
}
