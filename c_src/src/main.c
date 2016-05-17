#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "af_fourier.c"
#include <string.h>

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

double gradgrad(double x, double y, int K1, int K2)
{
	int m1 = K1/Nsize; int n1 = K1%Nsize;
	int m2 = K2/Nsize; int n2 = K2%Nsize;
	if(abs(m1-m2)>2 || abs(n1-n2)>2)
		return 0.;
	double dxdy = h/32., res;
	res =  0.25*(psi_b3(x+dxdy,y,K1)-psi_b3(x-dxdy,y,K1))*
			(psi_b3(x+dxdy,y,K2)-psi_b3(x-dxdy,y,K2))/dxdy +
		 0.25*(psi_b3(x,y+dxdy,K1)-psi_b3(x,y-dxdy,K1))*
			(psi_b3(x,y+dxdy,K2)-psi_b3(x,y-dxdy,K2))/dxdy;
	return res;
}

double integral_left(int K1, int K2)
{
	int m1 = K1/Nsize; int n1 = K1%Nsize;
	int m2 = K2/Nsize; int n2 = K2%Nsize;
	if(abs(m1-m2)>2 || abs(n1-n2)>2)
		return 0.;
	double dxdy = h/32., res=0.;
	int i,j;
	for(i=0; i<30; i++)
		for(j=0; j<30; j++)
		{
			res += dxdy*dxdy*gradgrad(h*(double)m1+dxdy*i, h*(double)n1+dxdy*j,K1,K2);
		}
	return res;
}

double integral_right(int K)
{
	double dxdy = h/32., res=0.;
	int i,j;
	int m1 = K/Nsize; int n1 = K%Nsize;
	for(i=0; i<30; i++)
		for(j=0; j<30; j++)
		{
			res += dxdy*dxdy*f(h*(double)m1+dxdy*i, h*(double)n1+dxdy*j)*
				psi_b3(h*(double)m1+dxdy*i, h*(double)n1+dxdy*j,K);
		}
	return res;
}

void matrix_solver()
{
	int i, j,k,l, NNsize = Nsize*Nsize;
	double pp =0.0;
	gsl_matrix * system 	= gsl_matrix_alloc (NNsize,NNsize);
	gsl_vector * coef		= gsl_vector_alloc (NNsize);
	gsl_vector * rightpart	= gsl_vector_alloc (NNsize);
	
	for(i=0; i<NNsize;i++)
	{
		gsl_vector_set(rightpart,i,integral_right(i));
		for(j=0; j<NNsize; j++)
		{
			gsl_matrix_set(system,i,j,integral_left(i,j));
		}
	}
	gsl_permutation * p = gsl_permutation_alloc (NNsize);
	gsl_linalg_LU_decomp (system, p, &i);
	gsl_linalg_LU_solve (system, p, rightpart, coef);
	
	FILE *op;
	op = fopen("plot.b3", "w");
	for(i=0; i<1024; i++)
	{
		for(j=0; j<1024; j++)
		{
			for(k=0;k<NNsize;k++)
				pp+=gsl_vector_get(coef,k)*psi_b3(i*(b-a)/1024.,j*(b-a)/1024.,k);
			fprintf(op,"%f %f %f\n",i*(b-a)/1024.,j*(b-a)/1024.,pp);
			pp = 0.;
		}
	}
	fclose(op);
}

/*			gsl_vector_set(b,0,0.0);*/
/*			gsl_matrix_set(sys, 0,0,f_B_3(-1.0));*/
/*			gsl_matrix_set(sys, 0,1,f_B_3(0.0));*/
/*			gsl_matrix_set(sys, 0,2,f_B_3(1.0));*/

void matrix_solver_v1()
{
	int i, j;
	int NNsize = Nsize*Nsize;
	
	gsl_matrix * system 	= gsl_matrix_alloc (NNsize,NNsize);
	gsl_vector * coef		= gsl_vector_alloc (NNsize);
	gsl_vector * rightpart	= gsl_vector_alloc (NNsize);
	
	
}





/*void solveB3(double A, double B, int n)*/
/*{*/
/*	int i,j;*/
/*	double outp, arg;*/
/*	const double step1=(B-A)/(n+1.0), step=(B-A)/(n-1.0), c0=f_dd_B_3(0.0)/step/step,c1=f_dd_B_3(1.0)/step/step ;*/
/*	*/
/*	gsl_matrix * sys = gsl_matrix_alloc (n+2, n+2);*/
/*	gsl_vector * x = gsl_vector_alloc (n+2);*/
/*	gsl_vector * b = gsl_vector_alloc (n+2);*/

/*	for(i=0;i<n+2;i++)*/
/*	{*/
/*		if(i==0)*/
/*		{*/
/*			gsl_vector_set(b,0,0.0);*/
/*			gsl_matrix_set(sys, 0,0,f_B_3(-1.0));*/
/*			gsl_matrix_set(sys, 0,1,f_B_3(0.0));*/
/*			gsl_matrix_set(sys, 0,2,f_B_3(1.0));*/

/*		}*/
/*		*/
/*		if(i!=0 && i!=n+1)*/
/*		{*/
/*			gsl_vector_set(b,i, f((double)(i-1)*step+A));*/
/*			gsl_matrix_set(sys,i,i-1,c1);*/
/*			gsl_matrix_set(sys,i,i,c0);*/
/*			gsl_matrix_set(sys,i,i+1,c1);*/
/*		}*/
/*		*/
/*		*/
/*		if(i==n+1)*/
/*		{*/
/*			gsl_vector_set(b,i,0.0);*/
/*			gsl_matrix_set(sys,n+1,n-1,f_B_3(-1.0));*/
/*			gsl_matrix_set(sys,n+1,n,f_B_3(0.0));*/
/*			gsl_matrix_set(sys,n+1,n+1,f_B_3(1.0));*/
/*		}	*/
/*	}*/
/*	//gsl_matrix_fprintf (stdout, sys, "%g");*/
/*	gsl_permutation * p = gsl_permutation_alloc (n+2);*/
/*	gsl_linalg_LU_decomp (sys, p, &i);*/
/*	gsl_linalg_LUsolve (sys, p, b, x);*/
/*	//gsl_vector_fprintf (stdout, x, "%g");*/
/*	*/
/*	FILE *op;*/
/*	op = fopen("./output/plot.b3", "w");*/
/*	for(arg = A; arg<=B; arg+=step)*/
/*	{*/
/*		outp = 0.0;*/
/*		for(i=0;i<n+2;i++)*/
/*			outp += gsl_vector_get(x,i)*f_B_3((arg-step*(double)(i-1))/step);*/
/*		fprintf(op,"%f %f\n",arg,outp);*/
/*	}*/
/*	fclose(op);*/
/*	op = fopen("./output/plot.b3.err", "w");*/
/*	for(arg = A; arg<=B; arg+=0.01)*/
/*	{*/
/*		outp = 0.0;*/
/*		for(i=0;i<n+2;i++)*/
/*			outp += gsl_vector_get(x,i)*f_B_3((arg-step*(double)(i-1))/step);*/
/*			*/
/*		fprintf(op,"%f %f\n",arg,fabs(outp-f_e(arg)));*/
/*	}*/
/*	fclose(op);*/

/*	*/
/*	op = fopen("./output/solution", "w");*/
/*	for(arg=A;arg<=B;arg+=0.01)*/
/*	{*/
/*		fprintf(op,"%f %f\n",arg,f_e(arg));*/
/*	}*/
/*	fclose(op);*/
/*	*/
/*}*/


int main (int arc, char** argv)
{
	matrix_solver();
	return 0;
}
