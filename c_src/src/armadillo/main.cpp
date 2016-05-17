#include <iostream>
#include <armadillo>
#include <cmath> 
#include "../B-splines.c"
using namespace std;
using namespace arma;


double a = -M_PI, b = M_PI, h = 2.*M_PI/(7.); 	//x&y a - min, b - max
int Nsize = 8;		//size by x and y,		//step - step by x&y
double graphics[1024][1024];



double f(double x, double y)
{
	return -2.*sin(x)*sin(y);
}

double omega(double x, double y)
{
	double ret = (M_PI*M_PI - x*x)*(M_PI*M_PI - y*y);
	if(ret >= 0.) return ret;
	return 0.;
}

double phi_b3(double x, double y)
{
	return f_B_3(x)*f_B_3(y);
}

double psi_b3(double x, double y, int K)
{
	int m = K/Nsize; int n = K%Nsize;
	double x0 = x - h*(double)m,
		 y0 = y - h*(double)n;
	if(fabs(x0/h)>1. || fabs(y0/h)>1. )
		return 0.;
	
	return omega(x,y)*phi_b3(2.*x0, 2.*y0);
}

double gradgrad()
{
	
}

double integral_left()
{
	
}

double integral_right()
{
	
}

void matrix_solver()
{
	
}


int main()
{
	mat A(5, 5, fill::randu);
	mat B = randu<mat>(4,5);
	cout << (A*inv(A)).t() << endl;
	return 0;
}
