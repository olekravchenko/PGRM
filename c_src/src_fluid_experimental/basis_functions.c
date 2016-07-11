double (*phi)(double, double, int);

double polynomial (double x, double y, int n)
{
    return pow(x,n%N)*pow(y,n/N);
}

double chebishov_1d(double x, int n)
{
    if(n == 0)
        return 1.;
    if(n == 1)
        return x;
    if(n == 2)
        return 2.*x*x-1.;
    if(n == 3)
        return 4.*x*x*x-3.*x;
    if(n == 4)
        return 8.*x*x*(x*x-1.)+1.;
    if(n == 5)
        return 4.*x*x*x*(4.*x*x-5.)+5.*x;
    if(n > 5)
        return 2.*x*chebishov_1d(x,n-1)-chebishov_1d(x,n-2);    
    return 0.;
}

double dchebishov_1d(double x, int n)
{
    if(n == 0)
        return 0.;
    if(n == 1)
        return 1.;
    if(n == 2)
        return 4.*x;
    if(n == 3)
        return 12.*x*x-3.;
    if(n == 4)
        return 32.*x*x*x-16.*x;
    if(n == 5)
        return 16.*5.*x*x*x*x-60.*x*x+5.;
   
    return 0.;
}

double chebishov_2d(double x, double y, int n)
{
    return chebishov_1d(x/(X1-X0),n%N)*chebishov_1d(y/(Y1-Y0),n/N);
}

double cubic_stepx, cubic_stepy;

double cubic_b_splines (double x, double y, int n)
{
    cubic_stepx = (X1-X0)/(double)(N-1);
    cubic_stepy = (Y1-Y0)/(double)(N-1);
    return f_B_3(0.5*cubic_stepx*(x-cubic_stepx*(double)(n%(N))))*
           f_B_3(0.5*cubic_stepy*(y-cubic_stepy*(double)(n/(N))));
}

double up_basis (double x, double y, int n)
{
    cubic_stepx = (X1-X0)/(double)(N-1);
    cubic_stepy = (Y1-Y0)/(double)(N-1);
    return f_up(0.5*cubic_stepx*(x-cubic_stepx*(double)(n%(N))))*
           f_up(0.5*cubic_stepy*(y-cubic_stepy*(double)(n/(N))));
}

//double fup_basis (double x, double y, int n)
//{
    //cubic_stepx = (X1-X0)/(double)(N-1);
    //cubic_stepy = (Y1-Y0)/(double)(N-1);
    //return f_cup(0.5*cubic_stepx*(x-cubic_stepx*(double)(n%(N))))*
           //f_cup(0.5*cubic_stepy*(y-cubic_stepy*(double)(n/(N))));
//}

void init_basis(int id)
{
    //	if(id == 1)
    phi = &cubic_b_splines;
    if(id == 2)
        phi = &polynomial;
    if(id == 3)
        phi = &chebishov_2d;
    if(id == 4)
	{	
		init_up();
		phi = &up_basis;
	}
}
