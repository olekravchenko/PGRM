#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>

double a = 0., b = M_PI, h = M_PI/4.; 	//x&y a - min, b - max
int Nsize = 25;		//size by x and y,		//step - step by x&y
double graphics[1024][1024];
double node[25][2];

double traps(double (*f)(double), double x0, double x1) 
//Trapezoidal rule integration for 1D
{
	double i;
	double step =(x1-x0)/50.;
	double ret = (*f)(x0);
	for(i=x0+step;i<x1; i+=step)
	{
		ret += 2.*(*f)(i);
	}
	ret *=0.5*step;
	return ret;
}

double integral2Dgrgr(double (*f)(double, double, int, int), double x0, double x1, double y0, double y1, int k1, int k2)
{
	double i,j, stepx=(x1-x0)/50., stepy=(y1-y0)/50.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{
			res += stepx*stepy*(*f)(i,j,k1,k2);
		}
	}
	return res;
}

double integral2Dfp(double (*f)(double, double, int), double x0, double x1, double y0, double y1, int k1)
{
	double i,j, stepx=(x1-x0)/50., stepy=(y1-y0)/50.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{
			res += stepx*stepy*(*f)(i,j,k1);
		}
	}
	return res;
}

double F(double x, double y)
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
	//int m = K/Nsize; int n = K%Nsize;
	double x0 = x - node[K][0],
		 y0 = y - node[K][1];
	//if(fabs(x0/h)>1. || fabs(y0/h)>1. )
	//	return 0.;
	
	return omega(x,y)*phi_b3(x0, y0);
}

double gradgrad(double x, double y, int K1, int K2)
{
	int m1 = K1/5; int n1 = K1%5;
	int m2 = K2/5; int n2 = K2%5;
	if(abs(m1-m2)>2 || abs(n1-n2)>2)
		return 0.;
	double dxdy = 0.02, res;
	res =  	0.25*(psi_b3(x+dxdy,y,K1)-psi_b3(x-dxdy,y,K1))*
			(psi_b3(x+dxdy,y,K2)-psi_b3(x-dxdy,y,K2))/dxdy +
		0.25*(psi_b3(x,y+dxdy,K1)-psi_b3(x,y-dxdy,K1))*
			(psi_b3(x,y+dxdy,K2)-psi_b3(x,y-dxdy,K2))/dxdy;
	return res;
}

double integral_left(int K1, int K2)
{
	double res;
	res = integral2Dgrgr(gradgrad,a,b,a,b,K1,K2);
	return res;
}

double fpsi(double x, double y, int k)
{
	return F(x,y)*psi_b3(x,y,k);
}


double integral_right(int K)
{
	return integral2Dfp(fpsi,a,b,a,b,K);
}

double reconstruct(gsl_vector *C, double x, double y)
{
	int i;
	double res=0.;
	
	for(i=0; i<Nsize; i++)
	{
		res += gsl_vector_get(C,i)*psi_b3(x,y,i);
	}
	
	return res;
}

void matrix_solver()
{
	int i, j,k,l;
	double pp =0.0;
	gsl_matrix * system 	= gsl_matrix_alloc (Nsize,Nsize);
	gsl_vector * coef	= gsl_vector_alloc (Nsize);
	gsl_vector * rightpart	= gsl_vector_alloc (Nsize);
	
	for(i=0; i<Nsize;i++)
	{
		//if(i<3 && j<3)
		gsl_vector_set(rightpart,i,-integral_right(i));
		for(j=0; j<Nsize; j++)
		{
			//if(i<3 && j<3)
			gsl_matrix_set(system,i,j,integral_left(i,j));
		}
	}
	FILE *op;
	op = fopen("matrix", "w");
	gsl_matrix_fprintf(op,system,"%g");
	fclose(op);
	
	gsl_permutation * p = gsl_permutation_alloc (Nsize);
	gsl_linalg_LU_decomp (system, p, &i);
	gsl_linalg_LU_solve (system, p, rightpart, coef);
	gsl_vector_fprintf(stdout,rightpart,"%g");
	printf("\n\n");
	gsl_vector_fprintf(stdout,coef,"%g");
	//FILE *op;
	op = fopen("plot.b3", "w");
	double X,Y;
	for(i=0; i<128; i++)
		for(j=0;j<128;j++)
		{
			X=a+(b-a)/64.*(double)(i-32.);
			Y=a+(b-a)/64.*(double)(j-32.);
			fprintf(op,"%f %f %f\n",X,Y,reconstruct(coef,X,Y));
		}
	fclose(op);
}


int main (int arc, char** argv)
{
	int i;	
	for(i = 0; i<25 ; i++)
	{
		node[i][0] = a + h*((double)(i%5));
		node[i][1] = a + h*((double)(i/5));
	//	printf("%f %f\n",node[i][0],node[i][1]);
	}
	matrix_solver();
	return 0;
}
