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
	double x0 = -1., x1 = 1., y0 = -1., y1 = 1.;
    double hx = (x1-x0)/512.,
           hy = (y1-y0)/512.,
           i,j;

    FILE * op;
    op = fopen("../../plot_data/plot_region", "w");
    for(i=x0; i<=x1; i+=hx)
        for(j=y0; j<=y1; j+=hy)
        {
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, (*f)(i,j));
        }
    fclose(op);
    i = system("../../bin/plotter_colorplot.py ../../plot_data/plot_region Numerical &");
}

int main()
{
	double omega_n (int n, double x, double y)
	{
		double p =2./3.;
		omega_primitive A = {.x0 = 0., .y0 = 0., .a = 1./3., .b = 1./3.};

		
		if(n == 0)
			return 1.;
		if(n == 1)
			return -rectangle(A, x, y, 0,0);
		
		
		
		return  R_and(-rectangle(A, x, y, 0,0),
				R_and( omega_n(n-1, 3.*(x - p), 3.*y),
				R_and( omega_n(n-1, 3.*(x - p), 3.*(y + p)),
				R_and( omega_n(n-1, 3.*x		, 3.*(y + p)),
				R_and( omega_n(n-1, 3.*(x + p), 3.*(y + p)),
				R_and( omega_n(n-1, 3.*(x + p), 3.*y		),
				R_and( omega_n(n-1, 3.*(x + p), 3.*(y - p)),
				R_and( omega_n(n-1, 3.*x		, 3.*(y - p)),
					   omega_n(n-1, 3.*(x - p), 3.*(y - p))))))))));
			
	}
	
	double omega_n2plot(double x, double y)
	{
		omega_primitive A = {.x0 = 0., .y0 = 0., .a = 1., .b = 1.};
		return fmax(0.,R_and(omega_n(4,x,y),rectangle(A, x, y, 0,0)));
	}
	plot_func(omega_n2plot);
	
	return 0;
}
