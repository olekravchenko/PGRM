
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
double integralLeft_trap(double (*f)(double,double,int,int), double x0, double x1, double y0, double y1, int k1, int k2)
{
	double i,j, stepx=(x1-x0)/100., stepy=(y1-y0)/100.;
	double res = 0.;
	
	for(i=x0+0.5*stepx; i<x1-stepy; i+= stepx)
	{
		res += 0.5*(*f)(i,y0,k1,k2);
		for(j=y0; j<y1; j+= stepy)
		{
			res += (*f)(i,j,k1,k2);
		}
		res += 0.5*(*f)(i,y1,k1,k2);
	}
	return res*stepx*stepy;
}
double integralRight_trap(double (*f)(double, double, int), double x0, double x1, double y0, double y1, int m)
{
	double i,j, stepx=(x1-x0)/100., stepy=(y1-y0)/100.;
	double res = 0.;
	
	for(i=x0+0.5*stepx; i<x1-stepy; i+= stepx)
	{
		res += 0.5*(*f)(i,y0,m);
		for(j=y0; j<y1; j+= stepy)
		{
			res += (*f)(i,j,m);
		}
		res += 0.5*(*f)(i,y1,m);
	}
	return res*stepx*stepy;
}
