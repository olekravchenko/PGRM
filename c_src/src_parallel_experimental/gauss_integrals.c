double nodes[4], weights[4];

void initGaussInt()
{
	//initializing node & weights
	nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
	nodes[0] = -nodes[3];
	nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
	nodes[1] = -nodes[2];
	
	weights[2] = 0.5+sqrt(30)/36;
	weights[3] = 0.5-sqrt(30)/36;
	weights[0] = weights[3];
	weights[1] = weights[2];
}

double SubIntegralLeft(double (*f)(double,double, int, int),
                       double x0,double x1,
                       double yc,
                       int k1,int k2)
{
	int i,j;//, N=32;
	double res = 0., step = (x1-x0)/intStep;
	//double x = 0.;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			//x = (double)(i-1)*step + a + 0.5*(nodes[j]+1.)*step;
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
	int i,j;//, N=32;
	double res = 0., step = (y1-y0)/intStep;
	//double x = 0.;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			//x = (double)(i-1)*step + a + 0.5*(nodes[j]+1.)*step;
			res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
		}
	}
	
	return 0.5*res*step;
}

double SubIntegralRight(double (*right_f)(double, double),
                        double (*sys_basis)(double, double, int),
                        double x0, double x1,
                        double yc, int k1)
{
	int i,j;//, N=32;
	double res = 0., step = (x1-x0)/intStep;
	//double x = 0.;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			//x = (double)(i-1)*step + a + 0.5*(nodes[j]+1.)*step;
			res += weights[j]
				*(*right_f)((double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step, yc)
				*(*sys_basis)((double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step, yc,k1);
		}
	}
	
	return 0.5*res*step;
}

//double (*f)(double,double, int),
double Simpson_right(double (*right_f)(double, double),
                     double (*sys_basis)(double, double, int),
                     double x0, double x1,
                     double y_c, int k1)
{
    double 	step = (x1-x0)/intStep;

    const int   q_of_layers = (int)((x1-x0)/step) + 1;// - (int)((x1-x0)/step)%2;
    int         k;

    double res =  (*right_f)(x0, y_c)*(*sys_basis)(x0, y_c, k1) +
                  (*right_f)((double)(q_of_layers-1)*step + x0, y_c)*(*sys_basis)((double)(q_of_layers-1)*step + x0, y_c, k1);
    //(*f)((double)(q_of_layers-1)*step + x0, y_c, k1);
    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*
          (*right_f)((double)k*step + x0, y_c)*(*sys_basis)((double)k*step + x0, y_c, k1);
        //(*f)((double)k*step + x0, y_c, k1);
    }
    return step/3.*res;
}


//simpson right integral
double integralRight_old(double (*right_f)(double, double),
                     double (*sys_basis)(double, double, int),
                     double x0, double x1,
                     double y0, double y1,
                     int k1)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
    double step = (y1-y0)/intStep;
    const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
    int         k;

    double res =  SubIntegralRight((*right_f),(*sys_basis),x0,x1,y0, k1) +
                  SubIntegralRight((*right_f),(*sys_basis),x0,x1,y0+(double)(q_of_layers-1)*step, k1);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*SubIntegralRight((*right_f),(*sys_basis),x0,x1,y0+(double)k*step, k1);
    }
    return step/3.*res;
}


double integralRight(double (*right_f)(double, double),
                     double (*sys_basis)(double, double, int),
                     double x0, double x1,
                     double y0, double y1,
                     int k1)
{
	int i,j;//, N=32;
	double res = 0., step = (y1-y0)/intStep;
	//double x = 0.;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			//x = (double)(i-1)*step + a + 0.5*(nodes[j]+1.)*step;
			res += weights[j]*SubIntegralRight((*right_f),(*sys_basis),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1);
		}
	}
	
	return 0.5*res*step;
}
