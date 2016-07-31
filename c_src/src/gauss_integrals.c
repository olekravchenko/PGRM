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

double SubIntegralRight(double (*right_f)(double, double),
                        double (*sys_basis)(double, double, int),
                        double x0, double x1,
                        double yc, int k1)
{
	int i,j;
	double res = 0., step = (x1-x0)/intStep;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			res += weights[j]
				*(*right_f)((double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step, yc)
				*(*sys_basis)((double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step, yc,k1);
		}
	}
	
	return 0.5*res*step;
}

double integralRight(double (*right_f)(double, double),
                     double (*sys_basis)(double, double, int),
                     double x0, double x1,
                     double y0, double y1,
                     int k1)
{
	int i,j;//, N=32;
	double res = 0., step = (y1-y0)/intStep;
	for (i = 1; i <= intStep; i++)
	{
		for (j = 0; j < 4; j++)
		{
			res += weights[j]*SubIntegralRight((*right_f),(*sys_basis),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1);
		}
	}
	
	return 0.5*res*step;
}
