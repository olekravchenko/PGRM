

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
