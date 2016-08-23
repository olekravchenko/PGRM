#include <stdbool.h>

double (*right_part_f)	(double, double);
double (*u_exact)		(double, double);
double (*f_boundary)	(double, double);
double (*omega)			(double, double);
double (*omega2)		(double, double);
double (*basis)			(double, double, int);
double (*phi)			(double, double, int);
double X0, X1, Y0, Y1;
int N;
//bool FiniteBasis = false;
double intStep, glob_delta, diff_step;

void init_eq(int);
