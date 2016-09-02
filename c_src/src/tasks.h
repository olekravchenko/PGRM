#include <stdbool.h>

/*
 * ToDo: 
 * - write normal comments, start writing documents for code
 * 
 */

double (*right_part_f)	(double, double);
//in case of Poisson's equation $\Delta u = f$, right_part_f is f itself

double (*f_boundary)	(double, double);
//function \varphi with Dirichlet values

double (*u_exact)		(double, double);
//exact solution of this problem

double (*omega)			(double, double);
//omega function for Dirichlet values (in mixed or first boundaty problem) or Neumann values (in second boundary problem)

double (*omega2)		(double, double);
//omega function for Neumann values (in mixed boundary problem only)

double (*structure)		(double, double, int);
//structure of solution. Usually is used one of function, defined in tasks.c

double (*phi)			(double, double, int);
//user-defined basis function with 'int' multi-index

double X0, X1, Y0, Y1;
//rectangle area of solution

int N;
double intStep, glob_delta, diff_step;

void init_eq(int);
