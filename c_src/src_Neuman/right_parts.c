double (*right_part_f)(double, double);
double (*u_exact)(double, double);
double X0, X1, Y0, Y1;
//int N;
#define Power pow
double f6(double x, double y)
{
	return (-16*x*(10*Power(x,8)*(35 + 9*Power(y,2)) + 
       Power(x,6)*(11575 - 5955*Power(y,2) + 
          405*Power(y,4) - 9*Power(y,6)) + 
       15*Power(x,4)*
        (-318655 + 35631*Power(y,2) - 
          1871*Power(y,4) + 23*Power(y,6)) + 
       5*Power(x,2)*
        (50504705 - 3918285*Power(y,2) + 
          99261*Power(y,4) - 887*Power(y,6) + 
          6*Power(y,8)) + 
       25*(-85942225 + 10099473*Power(y,2) - 
          188229*Power(y,4) - 
          1061*Power(y,6) + 42*Power(y,8))))/
   (675.*Power(1835 - 35*Power(y,2) + 
       Power(x,2)*(-35 + 3*Power(y,2)),3));
}

double f5(double x, double y) {return exp(x+y)*2. - 2.*exp(x+y)*((x-X0)*(x-X1)*(y-Y0)*(y-Y1)+(x-X0)*(x-X1)*(y-Y0)+(x-X0)*(x-X1)*(y-Y1)+(x-X0)*(x-X1)+(x-X0)*(y-Y0)*(y-Y1)+(x-X1)*(y-Y0)*(y-Y1)+(y-Y0)*(y-Y1));}
double u5(double x, double y) {return exp(x+y);}

double f4(double x, double y) {return 12.*(y*y*(x*x*x*x-625.) + x*x*(y*y*y*y-1.));}
double u4(double x, double y) {return (x*x*x*x-625.)*(y*y*y*y-1.);}

double f3(double x, double y) {return -2.*sin(x)*sin(y);}
double u3(double x, double y) {return sin(x)*sin(y);}

double f2(double x, double y) {return 2.*(x*x+y*y-2.);}
double u2(double x, double y) {return (x*x-1.)*(y*y-1.);}

double f1(double x, double y) {return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));}
double u1(double x, double y) {return (x*x*x*x-1.)*(y*y*y*y-1.);}

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
        X0 = 0.;
        X1 =  2.*M_PI;
        Y0 = 0.;
        Y1 =  2.*M_PI;
    }
    if(id == 4)
    {
        right_part_f = &f4;
        u_exact 	 = &u4;
        X0 = -5.;
        X1 =  5.;
        Y0 = -1.;
        Y1 =  1.;
    }
    if(id == 5)
    {
        right_part_f = &f5;
        u_exact 	 = &u5;
        X0 =  0.;
        X1 =  3.;
        Y0 =  0.;
        Y1 =  3.;
    }
    if(id == 6)
    {
        right_part_f = &f6;
        u_exact 	 = &u5;
        X0 =  -5.;
        X1 =  5.;
        Y0 =  -5.;
        Y1 =  5.;
    }
    
    
}