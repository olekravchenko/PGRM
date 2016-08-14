#include <stdio.h>
#include <stdlib.h>
int N=5;
double X0=-1.,X1=1.,Y0=-1.,Y1=1.;
#include "B-splines.c"
#include "af_poly.c"
#include "basis_functions.c"

int main(int argc, char** argv)
{
	init_basis(5);
    double hx = (X1-X0)/64.,
           hy = (Y1-Y0)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region", "w");
    for(i=X0; i<=X1; i+=hx)
        for(j=Y0; j<=Y1; j+=hy)
        {
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, phi(i,j,23));
		}
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region Numerical &");
	return 0;
}
