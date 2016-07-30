double integralLeft(double (*f)(double,double, int, int),
                        double x0,double x1,
                        double y0,double y1,
                        int k1,int k2)
{
    double stepy = (y1-y0)/intStep;
    double stepx = (x1-x0)/intStep;

    const int   q_y = (int)((y1-y0)/stepy) + 1;
	const int	q_x = (int)((x1-x0)/stepx) + 1; //-(int)((y1-y0)/step)%2;
    int         k,i;
	int 		kx = 1, ky = 1;
	double res = 0.;
    for(k = 0; k < q_y; k++)
    {
		if(k==0 || k==q_y-1)
			ky = 1;
		else if(k%2==1)
			ky = 4;
		else 
			ky = 2;
		for (i = 0; i < q_x; i++)
		{
			if(i==0 || i==q_x-1)
				kx = 1;
			else if(i%2==1)
				kx = 4;
			else
				kx = 2;
			res += ((double)ky*kx)*(*f)((double)i*stepx + x0, y0+(double)k*stepy, k1, k2);

				
		}
    }
    return stepx*stepy/9.*res;
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
double integralRight(double (*right_f)(double, double),
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
    
    double res =  Simpson_right((*right_f),(*sys_basis),x0,x1,y0, k1) +
		  Simpson_right((*right_f),(*sys_basis),x0,x1,y0+(double)(q_of_layers-1)*step, k1);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*Simpson_right((*right_f),(*sys_basis),x0,x1,y0+(double)k*step, k1);
    }
    return step/3.*res;
}
