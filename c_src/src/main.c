#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "B-splines.c"
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
#include "error_fs_new.c"



//#include "error_functions.c"
//ToDo: avoid usage of globally defined
//		re-write code in error_functions.c



double left_under_int_t(basis_args arguments, task Task)
{
    double 	x = arguments.x;
    double 	y = arguments.y;
    int 	m = arguments.m;
    int 	n = arguments.n;

    return  Task.Structure(x,y,m,Task)*(
            Task.Structure(x+diff_step,y,n,Task)+Task.Structure(x-diff_step,y,n,Task)+
            Task.Structure(x,y+diff_step,n,Task)+Task.Structure(x,y-diff_step,n,Task)
            -4.*Task.Structure(x,y,n,Task))*glob_delta*glob_delta;
}

double right_under_int_t(basis_args arguments, task Task)
{
    double 	x = arguments.x;
    double 	y = arguments.y;
    int 	m = arguments.m;

    return Task.right_part_f(x,y)*Task.Structure(x,y,m,Task);
}

void form_system_t (task *Task)
{
    int i, j;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;

    #pragma omp parallel for shared(Task,N) private(i,j) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral3(right_under_int_t,Task->area,args,2, Task));
        for(j = 0; j < N*N; j++)
        {
            args.n = j;
            gsl_matrix_set(Task->sys, i,j, gauss_integral3(left_under_int_t, Task->area, args,2, Task));
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

    #pragma omp parallel for shared(Task,N) private(i) firstprivate(args)
    for(i = 0; i < N*N; i++)
    {
        args.m = i;
        gsl_vector_set(Task->rightpart, i, gauss_integral3(right_under_int_t,Task->area,args,2, Task));
    }
}

void solve_task(task *Task)
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (Task->sys, p, &i);
    gsl_linalg_LU_solve (Task->sys, p, Task->rightpart, Task->solution);
}

void CFD_problem()
{
    double Reynolds_number = 10.;
	task stream_function, rotor_function;
    rect_area sol_area = {.x0 = X0, .x1 = X1, .y0 = Y0, .y1 = Y1};
    
    
    gsl_matrix *general_system = gsl_matrix_alloc (N*N,N*N); //single-use temporary storage for the matrix of the system

    //initial psi calculation
    tasks_constructor	(&stream_function,sol_area);
    form_system_t		(&stream_function);
    
    gsl_matrix_memcpy	(general_system, stream_function.sys);
    solve_task	(&stream_function);
    //plot_region_colorplot(stream_function.solution, stream_function.area);
    double psi(double x, double y)
    {
        return reconstruct_at_t(stream_function,x,y);
    }
    double laplacian(double (*f)(double, double),double x, double y)
    {
        return ((*f)(x+diff_step,y)+(*f)(x-diff_step,y)+
                (*f)(x,y+diff_step)+(*f)(x,y-diff_step)
                -4.*(*f)(x,y))*glob_delta*glob_delta;
    }
    double rotors_right_f(double x, double y)//ToDo: optimize this function, replace laplacian of laplacian with finite difference form
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
    double rotor_value	(double x, double y)
    {
        return reconstruct_at_t(rotor_function,x,y);
    }
    tasks_constructor	(&rotor_function,sol_area);

    rotor_function.f_boundary = 0;
    rotor_function.right_part_f = 0;
    rotor_function.Structure = 0;
    rotor_function.f_boundary = &rotor_boundary_f;
    rotor_function.right_part_f = &rotors_right_f;
    rotor_function.Structure = &Structure4;

    gsl_matrix_memcpy	(rotor_function.sys, general_system);
    form_right_part_t	(&rotor_function);
    solve_task	(&rotor_function);

    stream_function.right_part_f = 0;
    stream_function.right_part_f = &rotor_value;

    form_right_part_t	(&stream_function);
    gsl_matrix_memcpy	(rotor_function.sys, stream_function.sys);
    solve_task	(&stream_function);

    int i;
    for (i = 0; i < 120; i++)
    {
        form_right_part_t	(&rotor_function);
        solve_task	(&rotor_function);

        form_right_part_t	(&stream_function);
        gsl_matrix_memcpy	(rotor_function.sys, stream_function.sys);
        solve_task	(&stream_function);
		if(i%10 == 0)plot_lines_of_stream(stream_function.solution, stream_function.area);
    }

    plot_lines_of_stream(stream_function.solution, stream_function.area);
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
 */
{
    initGaussInt(16);

    omp_set_dynamic(1);
    omp_set_num_threads(8);
    int output_format = 0;

    if(argc>=6)
    {
        N 				= atoi(argv[1]);
        intStep 		= (double) atoi(argv[2]);
        init_eq		(atoi(argv[3]));
        init_basis	(atoi(argv[4]));
        output_format	= atoi(argv[5]);
    }
    else if(argc == 1)
    {

        N 			= 10;
        intStep 	= 1.;
        init_eq		(1);
        init_basis	(2);

        output_format	= 1000;
    }

    diff_step 	= pow(2.,-14);
    glob_delta 	= 1./diff_step;

	task function;// = &GetTask(1);
	//~ GetTask(1, &function);
    rect_area sol_area = {.x0 = X0, .x1 = X1, .y0 = Y0, .y1 = Y1};
    
    tasks_constructor	(&function,sol_area);
    form_system_t		(&function);
    
    //~ gsl_vector_fprintf(stdout,function.rightpart,"%f");
    
    solve_task	(&function);
    //~ multiplot(function.solution, function.area);
    basis_args zero_args = {0,0,0,0};
    
    multiplot_t(function);
    
	printf("%15.15e\n",errl2(zero_args,&function,2));
	
	
	//~ CFD_problem(); //to be used for testing solutions for Navier-Stokes equation in Stream function-Rotor form
	//task-id - 6 & 16
	
    return 0;
}
