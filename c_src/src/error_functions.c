//Depends on gauss_integrals initialization, if L_n error required




gsl_vector *solution_glob;



double error_L1_inner(double (*err_f)(double, double),
                      double x0, double x1,
                      double y_c)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
    double 	step = (x1-x0)/intStep;

    const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
    int         k;

    double res =  (*err_f)(x0, y_c) + (*err_f)((double)(q_of_layers-1)*step + x0, y_c);
    //(*f)((double)(q_of_layers-1)*step + x0, y_c, k1);
    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*(*err_f)((double)k*step + x0, y_c);
        //(*f)((double)k*step + x0, y_c, k1);
    }
    return step/3.*res;
}

double error_L1(double (*right_f)(double, double),
                double x0, double x1,
                double y0, double y1)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
    double step = (y1-y0)/intStep;
    const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
    int         k;

    double res =  error_L1_inner((*right_f),x0,x1,y0) +
                  error_L1_inner((*right_f),x0,x1,y0+(double)(q_of_layers-1)*step);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*error_L1_inner((*right_f),x0,x1,y0+(double)k*step);
    }
    return step/3.*res;
}

double error_L2_inner(double (*err_f)(double, double),
                      double x0, double x1,
                      double y_c)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
    double 	step = (x1-x0)/intStep;

    const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
    int         k;

    double res =  (*err_f)(x0, y_c)*(*err_f)(x0, y_c) +
                  (*err_f)((double)(q_of_layers-1)*step + x0, y_c)*(*err_f)((double)(q_of_layers-1)*step + x0, y_c),
                  temp;
    //(*f)((double)(q_of_layers-1)*step + x0, y_c, k1);
    for(k = 1; k < q_of_layers-1; k++)
    {
        temp = (*err_f)((double)k*step + x0, y_c);
        res += 2.*(double)(1+k%2)*temp*temp;
        //(*f)((double)k*step + x0, y_c, k1);
    }
    return step/3.*res;
}

double error_L2(double (*right_f)(double, double),
                double x0, double x1,
                double y0, double y1)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
    double step = (y1-y0)/intStep;
    const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
    int         k;

    double res =  error_L2_inner((*right_f),x0,x1,y0) +
                  error_L2_inner((*right_f),x0,x1,y0+(double)(q_of_layers-1)*step);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*error_L2_inner((*right_f),x0,x1,y0+(double)k*step);
    }
    return step/3.*res;
}

double error_f(double x, double y, gsl_vector *solution)
{
    return fabs(reconstruct_at(solution,x,y)-u_exact(x,y));
}

void errors_to_stdio(gsl_vector *solution,
                     double x1, double x2,
                     double y1, double y2)
{
    double  hx = (x2-x1)/64.,
            hy = (y2-y1)/64.,
            i,j,maxerr=0.,err;

    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            if((err = error_f(i,j))>maxerr) maxerr = err;
    printf("%d %g %f %f %e\n",N*N,intStep,
           maxerr,
           error_L1(error_f,x1,x2,y1,y2),
           error_L2(error_f,x1,x2,y1,y2));
}

double SubIntegralLeft(double (*f)(double,double, int, int),
                       double x0,double x1,
                       double yc,
                       int k1,int k2)
{
	int i,j;
	double res = 0., step = (x1-x0)/intStep;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			res += weights[j]*(*f)((double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step, yc, k1,k2);
		}
	}
	
	return 0.5*res*step;
}

double integralLeft(double (*f)(double,double, int, int),
                    double x0,double x1,
                    double y0,double y1,
                    int k1,int k2)
{
	int i,j;
	double res = 0., step = (y1-y0)/intStep;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
		}
	}
	
	return 0.5*res*step;
}
