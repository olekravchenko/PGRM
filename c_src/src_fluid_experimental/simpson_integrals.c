double Simpson_int(double (*f)(double,double, int, int),
                   double x0, double x1,
                   double y_c,
                   int k1, int k2)
{
    double	 step = (x1-x0)/intStep;
    
    const int   q_of_layers = (int)((x1-x0)/step) + 1;
    int         k;
    
    double res = (*f)(x0, y_c, k1, k2) +
                 (*f)((double)(q_of_layers-1)*step + x0, y_c, k1, k2);

    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*(*f)((double)k*step + x0, y_c, k1, k2);
    }
    return step/3.*res;
}
//Simpson_int_int
double integralLeft(double (*f)(double,double, int, int),
                    double x0,double x1,
                    double y0,double y1,
                    int k1,int k2)
{
    double step = (y1-y0)/intStep;
    
    const int   q_of_layers = (int)((y1-y0)/step) + 1; //-(int)((y1-y0)/step)%2;
    int         k;
    
    double res = Simpson_int((*f),x0,x1,y0, k1,k2) +
                 Simpson_int((*f),x0,x1,y0+(double)(q_of_layers-1)*step, k1,k2);
    
    for(k = 1; k < q_of_layers-1; k++)
    {
        res += 2.*(double)(1+k%2)*Simpson_int((*f),x0,x1,y0+(double)k*step, k1,k2);
    }
    return step/3.*res;
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
