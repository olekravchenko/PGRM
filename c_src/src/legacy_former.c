
double gauss_integral(	double (*f)(basis_args),
                        rect_area int_area,
                        basis_args args,
                        int dimension)
{
    double x0 = int_area.x0;
    double x1 = int_area.x1;
    double y0 = int_area.y0;
    double y1 = int_area.y1;

    int i,j, nodes_q = sizeof(nodes)/sizeof(double);
    double res = 0., stepx = (x1-x0)/intStep, stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < nodes_q; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
                res += weights[j]*gauss_integral((*f), int_area, temp_args, 1);
            }
        }

        return 0.5*res*stepy;
    }
    if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < nodes_q; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*(*f)(temp_args);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
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
