
__device__ double dev_f_B_3(double X)
{
	double absV=fabs(X);
	if(absV<2)
	{
		if(absV>=1)
			return 0.25*(2.0-absV)*(2.0-absV)*(2.0-absV); 
		else
			return 1.0 - 1.5*absV*absV*(1.0 - 0.5*absV);
	}
	return 0.0;
}

__host__ double host_f_B_3(double X)
{
	double absV=fabs(X);
	if(absV<2)
	{
		if(absV>=1)
		{
			return 0.25*(2.0-absV)*(2.0-absV)*(2.0-absV); 
		}
		else
			return 1.0 - 1.5*absV*absV*(1.0 - 0.5*absV);
	}
	return 0.0;
}


__device__ double (*phi)(double, double, int);

__device__ double phi_b3 (double x, double y, int n)
{
	//polynomial basis
//     return pow(x,n%N)*pow(y,n/N);
	
	double cubic_stepx = (X1-X0)/(double)(N-1);
	double cubic_stepy = (Y1-Y0)/(double)(N-1);
	
	return dev_f_B_3((N-1)/(X1-X0)*(x-X0-cubic_stepx*(double)(n%(N))))*
           dev_f_B_3((N-1)/(Y1-Y0)*(y-Y0-cubic_stepy*(double)(n/(N))));
}

__global__ void set_b3()
{
	phi = &phi_b3;
	
}

__device__ double omega(double x, double y) 
{
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}
__device__ double structure(double x, double y, int n)
{
	//structure for 1st boundary problem
    return phi(x,y,n)*omega(x,y);
	
	//structure for 2nd boundary problem
//     return 	phi(x,y,n)-omega(x,y)*
//             ((omega(x+diff_step,y)-omega(x-diff_step,y))*(phi(x+diff_step,y,n)-phi(x-diff_step,y,n))
//              +(omega(x,y+diff_step)-omega(x,y-diff_step))*(phi(x,y+diff_step,n)-phi(x,y-diff_step,n)))*
//             glob_delta*glob_delta*0.25;

}

__device__ double right_part_f(double x, double y)
{
    return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));
}



__host__ double Hphi (double x, double y, int n)
{
	//polynomial basis
//     return pow(x,n%N)*pow(y,n/N);

	double cubic_stepx = (X1-X0)/(double)(N-1);
	double cubic_stepy = (Y1-Y0)/(double)(N-1);
	
	return host_f_B_3((N-1)/(X1-X0)*(x-X0-cubic_stepx*(double)(n%(N))))*
           host_f_B_3((N-1)/(Y1-Y0)*(y-Y0-cubic_stepy*(double)(n/(N))));
	
}
__host__ double Homega(double x, double y) 
{
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}
__host__ double Hstructure(double x, double y, int n)
{
	//structure for 1st boundary problem
    return Hphi(x,y,n)*Homega(x,y);
	
	//structure for 2nd boundary problem
// 	return 	Hphi(x,y,n)-Homega(x,y)*
//             ((Homega(x+diff_step,y)-Homega(x-diff_step,y))*(Hphi(x+diff_step,y,n)-Hphi(x-diff_step,y,n))
//              +(Homega(x,y+diff_step)-Homega(x,y-diff_step))*(Hphi(x,y+diff_step,n)-Hphi(x,y-diff_step,n)))*
//             glob_delta*glob_delta*0.25;
}