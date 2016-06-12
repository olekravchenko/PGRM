
double integralLeft
     (double (*f)(double, double, int, int), 
      double x0, double x1, 
      double y0, double y1, 
      int k1, int k2)
// Returns value of left part integral: 
// I = \int_{x0}^{x1}\int_{y0}^{y1}\nabla\psi_{k1}\nabla\psi_{k2}dydx,
{
	double i,j, stepx=(x1-x0)/64., stepy=(y1-y0)/64.;
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

double integralRight
     (double (*f)(double, double, int), 
      double x0, double x1, 
      double y0, double y1, 
      int m)
// Returns value of integral in right part of equation:
// \int_{x0}^{x1}\int_{y0}^{y1} f(x,y)\psi_m(x,y)dydx
{
	double i,j, stepx=(x1-x0)/100., stepy=(y1-y0)/100.;
	double res = 0.;
	
	for(i=x0; i<x1; i+= stepx)
	{
		for(j=y0; j<y1; j+= stepy)
		{	
		      res += stepx*stepy*(*f)(i,j,m);
		}
	}
	return res;
}
