double (*right_part_f)(double, double);
double (*u_exact)(double, double);
double X0, X1, Y0, Y1;
//int N;


static double f3(double x, double y) {return -2.*sin(x)*sin(y);}
static double u3(double x, double y) {return sin(x)*sin(y);}

static double f2(double x, double y) {return 2.*(x*x+y*y-2.);}
static double u2(double x, double y) {return (x*x-1.)*(y*y-1.);}

static double f1(double x, double y) {return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));}
static double u1(double x, double y) {return (x*x*x*x-1.)*(y*y*y*y-1.);}

/*
static double (*f_s)(double, double);
static double (*u_s)(double, double);
*/

void init_eq(int id)
{
	if(id == 1)
	{
		right_part_f = &f1;
		u_exact 	 = &u1;
		X0 = -1.;
		X1 =  1.;
		Y0 = -1.;
		Y1 =  1.;
	}
	
	if(id == 2)
	{
		right_part_f = &f2;
		u_exact 	 = &u2;
		X0 = -1.;
		X1 =  1.;
		Y0 = -1.;
		Y1 =  1.;
	}
	
	if(id == 3)
	{
		right_part_f = &f3;
		u_exact 	 = &u3;
		X0 = -M_PI;
		X1 =  M_PI;
		Y0 = -M_PI;
		Y1 =  M_PI;
	}
}
