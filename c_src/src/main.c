#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_poly.c"
#include <string.h>
#include "smplDynArray.c"
#include "tasks.h"

#include "types_and_structs.c"

//Definitions moved to tasks.h
//ToDo: reduce usage of global variables
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
{
    int i, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
    
#pragma omp parallel for shared(system, RightPart,N) private(i,j) firstprivate(args)
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


double left_under_int_t(basis_args arguments, task Task)
{
	double 	x = arguments.x;
	double 	y = arguments.y;
	int 	m = arguments.m;
	int 	n = arguments.n;
	
    return  Task.structure(x,y,m)*(
			Task.structure(x+diff_step,y,n)+Task.structure(x-diff_step,y,n)+
			Task.structure(x,y+diff_step,n)+Task.structure(x,y-diff_step,n)
			-4.*Task.structure(x,y,n))*glob_delta*glob_delta;
}

double right_under_int_t(basis_args arguments, task Task)
{
	double 	x = arguments.x;
	double 	y = arguments.y;
	int 	m = arguments.m;
	
	return Task.right_part_f(x,y)*Task.structure(x,y,m);
}

void form_system_t (task *Task)
{	
    int i, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
    
//#pragma omp parallel for shared(system, RightPart,N) private(i,j) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral2(right_under_int_t,Task->area,args,2, Task));
        for(j = 0; j < N*N; j++)
        {
            args.n = j;
            gsl_matrix_set(Task->sys, i,j, gauss_integral2(left_under_int_t, Task->area, args,2, Task));
        }
    }
}

void form_right_part_t (task *Task)
{	
    int i;//, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
    
//#pragma omp parallel for shared(system, RightPart,N) private(i,j) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral2(right_under_int_t,Task->area,args,2, Task));
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
void solve_matrix_eq_t(task *Task)
//Solve SLE Ax=b, where A = system, b = RightPart, x = solution
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (Task->sys, p, &i);
    gsl_linalg_LU_solve (Task->sys, p, Task->rightpart, Task->solution);
}

double F3(double x, double y) {return -2.*sin(x)*sin(y);}
double rectangle(double x, double y){return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);}
double bf3(double x, double y){return 0.;}
    task stream_function, rotor_function;
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
 * 	-finish rewriting code for solving systems of PDE
 * 	-finish rewriting code with methods of description of OOP
 */
{
    initGaussInt();

    omp_set_dynamic(1);
    omp_set_num_threads(16);
    int output_format = 0;

    if(argc>=6)
    {
        N 				= atoi(argv[1]);
        intStep 		= (double) atoi(argv[2]);
        init_eq(atoi(argv[3]));
        init_basis(atoi(argv[4]));
        output_format	= atoi(argv[5]);
    }
    else if(argc == 1)
    {
        N 			= 8;
        intStep 	= 4.;
        init_eq(6);
        init_basis(3);
        output_format	= 1000;
    }

    diff_step 	= pow(2.,-9);
    glob_delta 	= 1./diff_step;
	
	
	double Reynolds_number = 10.;
	
    rect_area sol_area = {.x0 = X0, .x1 = X1, .y0 = Y0, .y1 = Y1};
    gsl_matrix *general_system = gsl_matrix_alloc (N*N,N*N); //temporary storage for the matrix of the system

    //initial psi calculation
    tasks_constructor	(&stream_function,sol_area);

    form_system_t		(&stream_function);
    gsl_matrix_memcpy	(general_system, stream_function.sys);
    solve_matrix_eq_t	(&stream_function);
    plot_by_argument	(stream_function.solution, output_format, stream_function.area);
    /*


        //init_eq(6);
        //init_basis(5);

    	//sol_area.x0 = X0;
    	//sol_area.x1 = X1;
    	//sol_area.y0 = Y0;
    	//sol_area.y1 = Y1;
    */
    //argc = system("sleep 1"); //using argc to not define new variable

    double psi(double x, double y)
    {
		//printf("%f %f\n",x ,y);
        return reconstruct_at_t(stream_function,x,y);
    }
    double laplacian(double (*f)(double, double),double x, double y)
    {
		//printf("%f %f\n",x ,y);
        return ((*f)(x+diff_step,y)+(*f)(x-diff_step,y)+
                (*f)(x,y+diff_step)+(*f)(x,y-diff_step)
                -4.*(*f)(x,y))*glob_delta*glob_delta;
    }
	double rotors_right_f(double x, double y)//ToDo: optimize this function
	{
		double	l_psi_px = laplacian(psi,x+diff_step,y);
		double	l_psi_py = laplacian(psi,x-diff_step,y);
		double	l_psi_mx = laplacian(psi,x,y+diff_step);
		double	l_psi_my = laplacian(psi,x,y-diff_step);
		double 	psi_px = psi(x+diff_step,y);
		double	psi_py = psi(x,y+diff_step);
		double	psi_mx = psi(x+diff_step,y);
		double	psi_my = psi(x,y+diff_step);
				
		return Reynolds_number*0.25*(
				((psi_py-psi_my)*(l_psi_px-l_psi_mx) -
				 (psi_px-psi_mx)*(l_psi_py-l_psi_my)) +
				(l_psi_px+l_psi_mx+l_psi_py+l_psi_my-4.*laplacian(psi,x,y)))*glob_delta*glob_delta;
	}
	double rotor_boundary_f(double x, double y) //ToDo: reduce one of this functions
	{
		return -laplacian(psi,x,y);
	}

    tasks_constructor	(&rotor_function,sol_area);
    gsl_matrix_memcpy	(rotor_function.sys, general_system);

	rotor_function.f_boundary = 0;
    rotor_function.right_part_f = 0;
	//f_boundary = 0;
    
	rotor_function.f_boundary = &rotor_boundary_f;
    rotor_function.right_part_f = &rotors_right_f;
	//f_boundary = &rotor_boundary_f;

    
	printf("1\n");
    form_right_part_t	(&rotor_function);
    printf("2\n");
    solve_matrix_eq_t	(&rotor_function);
    //plot_by_argument	(rotor_function.solution, output_format, rotor_function.area);

	double rotor_value	(double x, double y)
	{
        return -reconstruct_at_t(rotor_function,x,y);
	}
	stream_function.right_part_f = 0;
	stream_function.right_part_f = &rotor_value;
	
    form_right_part_t	(&stream_function);
    gsl_matrix_memcpy	(general_system, stream_function.sys);
    solve_matrix_eq_t	(&stream_function);
    plot_by_argument	(stream_function.solution, output_format, stream_function.area);
    return 0;
}
