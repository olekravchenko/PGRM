double (*right_part_f)	(double, double);
double (*u_exact)		(double, double);
double (*f_boundary)	(double, double);
double (*omega)			(double, double);
double (*omega2)		(double, double);
double (*basis)			(double, double, int);
double (*phi)			(double, double, int);
double X0, X1, Y0, Y1;

//Dirichlet problem structure
double basis1(double x, double y, int n)
{
    return phi(x,y,n)*omega(x,y);
}
//Neumann problem structure
double basis2(double x, double y, int n)
{
	return 	phi(x,y,n)-omega(x,y)*
			((omega(x+diff_step,y)-omega(x-diff_step,y))*(phi(x+diff_step,y,n)-phi(x-diff_step,y,n))
			+(omega(x,y+diff_step)-omega(x,y-diff_step))*(phi(x,y+diff_step,n)-phi(x,y-diff_step,n)))*
			glob_delta*glob_delta*0.25;
}
//Mixed boundary problem
double basisM(double x, double y, int n)
{
	return 	basis1(x,y,n)-omega(x,y)*omega2(x,y)/(omega(x,y)+omega2(x,y))*
			((omega2(x+diff_step,y)-omega2(x-diff_step,y))*(basis1(x+diff_step,y,n)-basis1(x-diff_step,y,n))
			+(omega2(x,y+diff_step)-omega2(x,y-diff_step))*(basis1(x,y+diff_step,n)-basis1(x,y-diff_step,n)))*
			glob_delta*glob_delta*0.25;
}
#define Power pow
#define Sqrt sqrt





double o10_2(double x, double y)
{
	return 0.4 + 2.5*x*x + 2.*y*y - 
   2.5*Sqrt(0.0416 + x*x*x*x - 0.16*y*y + 1.04*y*y*y*y + 
      x*x*(-0.08 + 2.*y*y));
}
double bf10(double x, double y)
{
	return x - (6.25*x*(-1. + x*x)*
      (0.19999999999999996 + 1.25*x*x + 1.*y*y - 
        1.25*Sqrt(0.0416 + x*x*x*x - 0.16*y*y + 1.04*y*y*y*y + 
           x*x*(-0.08 + 2.*y*y)))*
      (-0.04 + 1.*x*x + 1.*y*y - 
        1.*Sqrt(0.0416 + x*x*x*x - 0.16*y*y + 1.04*y*y*y*y + 
           x*x*(-0.08 + 2.*y*y))))/
    (Sqrt(0.25*Power(-1. + y*y,2) + 
        (25*Power(-0.04 + x*x + y*y,2))/4.)*
      (0.45 + 1.*x*x + 1.*y*y - 
        1.25*Sqrt(0.0416 + x*x*x*x - 0.16*y*y + 1.04*y*y*y*y + 
           x*x*(-0.08 + 2.*y*y))));
}


double o9(double x, double y)
{
	return 0.5*(1.-x*x);
}
double o9_2(double x, double y)
{
	return 0.5*(1.-y*y);
}
double bf9(double x, double y)
{
	return x;
}


double f_num(double x, double y)
{
	return -(f_boundary(x+diff_step,y)+f_boundary(x-diff_step,y)+
			f_boundary(x,y+diff_step)+f_boundary(x,y-diff_step)-
			4.*f_boundary(x,y))*glob_delta*glob_delta;
}

double f7(double x, double y)
{
	return (x*(-8*x*x*(-1 + y*y)*(-2 + x*x + y*y)*
        (1 + (-1 + x*x)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       4*(-1 + y*y)*Power(-2 + x*x + y*y,2)*
        (1 + (-1 + x*x)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       (-1 + y*y)*Power(-2 + x*x + y*y,2)*
        (2 - (4*Power(-x + Power(x,3),2))/
           Power(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y,1.5) + 
          (-2 + 6*x*x)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) - 
       8*y*y*(-1 + y*y)*(-2 + x*x + y*y)*
        (1 + (-1 + y*y)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       8*y*y*Power(-2 + x*x + y*y,2)*
        (1 + (-1 + y*y)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       (-1 + y*y)*Power(-2 + x*x + y*y,2)*
        (2 - (4*Power(-y + Power(y,3),2))/
           Power(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y,1.5) + 
          (-2 + 6*y*y)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       8*x*x*(-1 + y*y)*
        (-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       8*y*y*(-1 + y*y)*
        (-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) - 
       8*y*y*(-2 + x*x + y*y)*
        (-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) - 
       8*(-1 + y*y)*(-2 + x*x + y*y)*
        (-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)) + 
       2*Power(-2 + x*x + y*y,2)*
        (-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y))))/
   (2.*Power(-2 + x*x + y*y,3));
}
double bf7(double x, double y)
{
	return -(x*(-1 + y*y)*(-2 + x*x + y*y + 
        Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)))/
   (2.*(-2 + x*x + y*y));
}
double o7(double x, double y)
{
	return (2 - x*x - y*y - Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
       y*y*y*y))/2.;
}

double f6(double x, double y)
{
	return -(225.*y*(405000.*Power((-2.*x*(-1. + x*x))/
           Sqrt(2. - 2.*x*x + x*x*x*x - 2.*y*y + y*y*y*y) + 
          (-8.*x + 450.*x*x*x)/
           Sqrt(32. - 1800.*x*x + 50625.*x*x*x*x - 1800.*y*y + 
             50625.*y*y*y*y),2)*(-0.035555555555555556 + x*x + y*y + 
          Sqrt(32. - 1800.*x*x + 50625.*x*x*x*x - 1800.*y*y + 
             50625.*y*y*y*y)/225.) + 
       405000.*Power((-2.*y*(-1. + y*y))/
           Sqrt(2. - 2.*x*x + x*x*x*x - 2.*y*y + y*y*y*y) + 
          (-8.*y + 450.*y*y*y)/
           Sqrt(32. - 1800.*x*x + 50625.*x*x*x*x - 1800.*y*y + 
             50625.*y*y*y*y),2)*(-0.035555555555555556 + x*x + y*y + 
          Sqrt(32. - 1800.*x*x + 50625.*x*x*x*x - 1800.*y*y + 
             50625.*y*y*y*y)/225.) - 
       225.*((16.*Power(-x + x*x*x,2))/
           Power(2. - 2.*x*x + x*x*x*x - 2.*y*y + y*y*y*y,1.5) + 
          (8. - 24.*x*x)/
           Sqrt(2. - 2.*x*x + x*x*x*x - 2.*y*y + y*y*y*y) - 
          (3600.*Power(-4.*x + 225.*x*x*x,2))/
           Power(32. - 1800.*x*x + 50625.*x*x*x*x - 1800.*y*y + 
             50625.*y*y*y*y,1.5) + 
          (8.*(-4. + 675.*x*x))/
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y))*(-0.035555555555555556 + x*x + y*y + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y)/225.)*
        (442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)
          ) - 225*((16*Power(-y + y*y*y,2))/
           Power(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y,1.5) + 
          (8 - 24*y*y)/
           Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) - 
          (3600*Power(-4*y + 225*y*y*y,2))/
           Power(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y,1.5) + 
          (8*(-4 + 675*y*y))/
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y))*(-0.035555555555555556 + x*x + y*y + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y)/225.)*
        (442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)
          ) + 4*(2 - (900*Power(-4*x + 225*x*x*x,2))/
           Power(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y,1.5) + 
          (2*(-4 + 675*x*x))/
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y))*Power(442 - 
          225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)
          ,2) + 4*(2 - (900*Power(-4*y + 225*y*y*y,2))/
           Power(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y,1.5) + 
          (2*(-4 + 675*y*y))/
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y))*Power(442 - 
          225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)
          ,2) + (16*(-4 + 225*y*y + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*Power(442 - 
            225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y),2))/
        Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)\
        + (7200*x*x*(-4 + 225*x*x + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(442 - 
            225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(4*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
               y*y*y*y) - Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
              1800*y*y + 50625*y*y*y*y) + 
            x*x*(-225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
                  y*y*y*y) + Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
                 1800*y*y + 50625*y*y*y*y))))/
        (Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)*
          (32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)) + 
       (7200*y*y*(-4 + 225*y*y + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(442 - 
            225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(4*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
               y*y*y*y) - Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
              1800*y*y + 50625*y*y*y*y) + 
            y*y*(-225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
                  y*y*y*y) + Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
                 1800*y*y + 50625*y*y*y*y))))/
        (Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)*
          (32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)) + 
       (16*(-8 + 225*x*x + 225*y*y + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(442 - 
            225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
            Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
              50625*y*y*y*y))*(4*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
               y*y*y*y) - Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
              1800*y*y + 50625*y*y*y*y) + 
            y*y*(-225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
                  y*y*y*y) + Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 
                 1800*y*y + 50625*y*y*y*y))))/
        (Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y)*
          Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)
          )))/(4.*Power(442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + 
          y*y*y*y) + Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
         50625*y*y*y*y),3));
}
double bf6(double x, double y)
{
	return (y*(-8 + 225*x*x + 225*y*y + 
       Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y)))/
   (442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
     Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y));
}
double o6(double x, double y)
{
	return (442 - 225*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
     Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 50625*y*y*y*y) - 
     225*Sqrt(Power(-2 + x*x + y*y + 
          Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y),2) + 
        Power(-8 + 225*x*x + 225*y*y + 
           Sqrt(32 - 1800*x*x + 50625*x*x*x*x - 1800*y*y + 
             50625*y*y*y*y),2)/50625.))/225.;
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

double laplace_f(double x, double y){return 0.;}
double omega_rectangle(double x, double y){return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);}
void init_eq(int id)
{
	basis = &basis1;
    if(id == 1)
    {
        right_part_f = &f1;
        u_exact 	 = &u1;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        omega2		 = &omega_rectangle;
        X0 = -1.;
        X1 =  1.;
        Y0 = -1.;
        Y1 =  1.;
    }

    if(id == 2)
    {
        right_part_f = &f2;
        u_exact 	 = &u2;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        X0 = -1.;
        X1 =  1.;
        Y0 = -1.;
        Y1 =  1.;
    }

    if(id == 3)
    {
        right_part_f = &f3;
        u_exact 	 = &u3;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        omega2		 = &omega_rectangle;
        X0 = 0.;
        X1 =  2.*M_PI;
        Y0 = 0.;
        Y1 =  2.*M_PI;
    }
    if(id == 4)
    {
        right_part_f = &f4;
        u_exact 	 = &u4;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        X0 = -5.;
        X1 =  5.;
        Y0 = -1.;
        Y1 =  1.;
    }
    if(id == 5)
    {
        right_part_f = &f5;
        u_exact 	 = &u5;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        omega2		 = &omega_rectangle;
        X0 =  0.;
        X1 =  3.;
        Y0 =  0.;
        Y1 =  3.;
    }
    if(id == 6)
    {
        right_part_f = &f6;
        u_exact 	 = &u5;
        f_boundary	 = &bf6;
        omega		 = &o6;
        omega2		 = &o6;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }    
    if(id == 7)
    {
		basis		 = &basis2;
        right_part_f = &f7;
        u_exact 	 = &u5;
        f_boundary	 = &bf7;
        omega		 = &o7;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }    
    if(id == 8)
    {
		basis		 = &basis2;
        right_part_f = &f_num;
        u_exact 	 = &u5;
        f_boundary	 = &bf7;
        omega		 = &o7;
        omega2		 = &o7;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }    
	if(id == 9)
	{
		basis		 = &basisM;
        right_part_f = &f_num;
        u_exact 	 = &u5;
        f_boundary	 = &bf9;
        omega		 = &o9;
        omega2		 = &o9_2;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;	
    }
	if(id == 10)
	{
		basis		 = &basisM;
        right_part_f = &f_num;
        u_exact 	 = &u5;
        f_boundary	 = &bf10;
        omega		 = &o9;
        omega2		 = &o10_2;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
}