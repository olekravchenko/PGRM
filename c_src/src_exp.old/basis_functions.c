#include <math.h>
int N;
double X0, X1, Y0, Y1;
double (*phi)(double, double, int);

static double polynomial (double x, double y, int n)
{
	return pow(x,n%N)*pow(y,n/N);
}

static double cubic_b_splines (double x, double y, int n)
{
	return f_B_3(0.5*(x-(double)(n%N))/(X1-X0))*f_B_3(0.5*(y-(double)(n/N))/(Y1-Y0));
}

void init_basis(int id)
{
//	if(id == 1)
		phi = &cubic_b_splines;
	if(id == 2)
		phi = &polynomial;
}
