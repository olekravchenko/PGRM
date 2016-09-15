
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






















double left_under_int(double x, double y, int m, int n)
{
 // \phi_m \Delta \phi_n 
    return  structure(x,y,m)*(
			structure(x+diff_step,y,n)+structure(x-diff_step,y,n)+
			structure(x,y+diff_step,n)+structure(x,y-diff_step,n)
			-4.*structure(x,y,n))*glob_delta*glob_delta;
}

/*
 * 
 * Place for Including of parallel former, if required
 * 
 * //ToDo: 
 * 	-re-write, using new new definition of left and right under integral functions
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
        gsl_vector_set(RightPart, i, integralRight(right_part_f,structure,x1,x2,y1,y2,i));
        for(j = 0; j < N*N; j++)
        {
            gsl_matrix_set(system, i,j, integralLeft(left_under_int,x1,x2,y1,y2,i,j));
        }
    }
}
