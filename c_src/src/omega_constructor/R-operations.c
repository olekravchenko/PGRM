#include "main_omega.c"




//~ double R_and(double x, double y)
//~ {
	//~ return x+y-sqrt(x*x+y*y);
//~ }
//~ double R_or(double x, double y)
//~ {
	//~ return x+y+sqrt(x*x+y*y);
//~ }


void plot_func(double (*f)(double, double))
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
	double x0 = -2., x1 = 2., y0 = -2., y1 = 2.;
    double hx = (x1-x0)/128.,
           hy = (y1-y0)/128.,
           i,j;

    FILE * op;
    op = fopen("../../plot_data/plot_region", "w");
    for(i=x0; i<=x1; i+=hx)
        for(j=y0; j<=y1; j+=hy)
        {
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, (*f)(i,j));
        }
    fclose(op);
    i = system("../../bin/plotter.py ../../plot_data/plot_region Numerical &");
}

int main()
{
	omega_primitive A = {.x0 = -1., .y0 = -1., .a = 1., .b = 1.};
	double omega_new(double x, double y)
	{
		return line_segment(A,x,y);
	}
	plot_func(omega_new);
	
	return 0;
}
