#include "main_omega.c"




double R_and(double x, double y)
{
	return x+y-sqrt(x*x+y*y);
}
double R_or(double x, double y)
{
	return x+y+sqrt(x*x+y*y);
}


void plot_func(double (*f)(double, double))
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
	double x0 = -1., x1 = 1., y0 = -1., y1 = 1.;
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
	omega_primitive A = {.x0 = 0., .y0 = 0., .a = 1./5., .b = 2./5};
	double figure1(double x, double y)
	{
		return band_x(A, x,y);
	}
	double figure2(double x, double y)
	{
		return band_y(A, x,y);
	}
	
	double omega_new(double x, double y)
	{
		return fmax(0.,R_and(-band_x(A, x,y), band_y(A, x,y)));
	}
	plot_func(omega_new);
	
	return 0;
}
