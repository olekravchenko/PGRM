#include <math.h>
#include "tasks.h"
#include "omega_constructor/main_omega.c"
//Dirichlet problem structure
double structure1(double x, double y, int n)
{
    return phi(x,y,n)*omega(x,y);//*omega_at;
}
//Rigid clamping :-)
double structure4(double x, double y, int n)
{
    double omega_at = omega(x,y);
    return phi(x,y,n)*omega_at*omega_at;
}
//Neumann problem structure
double structure2(double x, double y, int n)
{
    return 	phi(x,y,n)-omega(x,y)*
            ((omega(x+diff_step,y)-omega(x-diff_step,y))*(phi(x+diff_step,y,n)-phi(x-diff_step,y,n))
             +(omega(x,y+diff_step)-omega(x,y-diff_step))*(phi(x,y+diff_step,n)-phi(x,y-diff_step,n)))*
            glob_delta*glob_delta*0.25;
}
//Mixed boundary problem
double structureM(double x, double y, int n)
{
    return 	structure1(x,y,n)-omega(x,y)*omega2(x,y)/(omega(x,y)+omega2(x,y))*
            ((omega2(x+diff_step,y)-omega2(x-diff_step,y))*(structure1(x+diff_step,y,n)-structure1(x-diff_step,y,n))
             +(omega2(x,y+diff_step)-omega2(x,y-diff_step))*(structure1(x,y+diff_step,n)-structure1(x,y-diff_step,n)))*
            glob_delta*glob_delta*0.25;
}

//Newton/Robin/3rd order		ToDo: fix it and test it
double structure3(double x, double y, int n)
{
	double h=0.;
    return 	phi(x,y,n)-omega(x,y)*(
            (-h+omega(x+diff_step,y)-omega(x-diff_step,y))*(phi(x+diff_step,y,n)-phi(x-diff_step,y,n))
             +(omega(x,y+diff_step)-omega(x,y-diff_step))*(phi(x,y+diff_step,n)-phi(x,y-diff_step,n)))*
            glob_delta*glob_delta*0.25;
}


#define Power pow
#define Sqrt sqrt

double bf17(double x, double y)
{
	omega_primitive external = {.x0=0., .y0 = 0., .a=X1, .b=Y1};
	omega_primitive inner	 = {.x0=x_center, .y0 = y_center, .a=X1/5., .b=Y1/5.};

	double	inner_v = rectangle(inner, x,y,0,0), 
			extern_v = rectangle(external, x,y,0,0);
	return (/*0./extern_v*/ - 1./inner_v)/(-1./inner_v + 1./extern_v);
}
double o17(double x, double y)
{
	omega_primitive external = {.x0=0., .y0 = 0., .a=X1, .b=Y1};
	omega_primitive inner	 = {.x0=x_center, .y0 = y_center, .a=X1/5., .b=Y1/5.};
	return R_and(rectangle(external, x,y,0,0),-rectangle(inner, x,y,0,0));
}



double bf16(double x, double y)
{
	omega_primitive external = {.x0=0., .y0 = 0., .a=X1, .b=Y1};
	omega_primitive inner	 = {.x0=0., .y0 = 0., .a=X1/5., .b=Y1/5.};

	double	inner_v = rectangle(inner, x,y,0,0), 
			extern_v = rectangle(external, x,y,0,0);
	return (inner_v*y/*(1.-y*y/3.)*/-extern_v*0.)/(inner_v - extern_v);
}
double o16(double x, double y)
{
	omega_primitive external = {.x0=0., .y0 = 0., .a=X1, .b=Y1};
	omega_primitive inner	 = {.x0=0., .y0 = 0., .a=X1/5., .b=Y1/5.};
	return R_and(rectangle(external, x,y,0,0),-rectangle(inner, x,y,0,0));
}


double bf15(double x, double y)
{
	
	return ((-0.75 - 5.333333333333333*Power(-0.75 + y,3) + y)/(1.5 + x + Sqrt((1+x)*(1+x) + (y-0.5)*(y-0.5)) - 1.*y) + 
     (0.75 + y - 5.333333333333333*Power(0.75 + y,3))/(1.5 - 1.*x + y + Sqrt(1.25 - 2*x + x*x + 1.*y + y*y)) - 
     0.16666666666666666/(-1.5 + x + 2.*y + Sqrt(2 - 2*x + x*x - 2*y + y*y) - 
        1.*Sqrt((y+0.5)*(y+0.5) + Power(-2 + x + y + Sqrt(2 - 2*x + x*x - 2*y + y*y),2))) - 
     1/(6.*(1.5 + x + 2*y - Sqrt((1+x)*(1+x) + (1+y)*(1+y)) + 
          Sqrt((y-0.5)*(y-0.5) + Power(2 + x + y - Sqrt((1+x)*(1+x) + (1+y)*(1+y)),2)))))/
   (1/(1.5 + x + Sqrt((1+x)*(1+x) + (y-0.5)*(y-0.5)) - 1.*y) + 
     1/(1.5 - 1.*x + y + Sqrt(1.25 - 2*x + x*x + 1.*y + y*y)) - 
     1./(-1.5 + x + 2.*y + Sqrt(2 - 2*x + x*x - 2*y + y*y) - 
        1.*Sqrt((y+0.5)*(y+0.5) + Power(-2 + x + y + Sqrt(2 - 2*x + x*x - 2*y + y*y),2))) + 
     1/(1.5 + x + 2*y - Sqrt((1+x)*(1+x) + (1+y)*(1+y)) + 
        Sqrt((y-0.5)*(y-0.5) + Power(2 + x + y - Sqrt((1+x)*(1+x) + (1+y)*(1+y)),2))));
}

double o14(double x, double y)
{
	omega_primitive external = {.x0=0., .y0 = 0., .a=X1, .b=Y1};
	omega_primitive inner	 = {.x0=0., .y0 = 0., .a=X1/25., .b=Y1/25.};

    return R_and(rectangle(external, x,y,0,0),-rectangle(inner, x,y,0,0));
}
double bf14(double x, double y)
{
    return -10.*((y*(-3 + x*x + y*y + 
         Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y))*
       (-2 + 625*x*x + 625*y*y + 
         Sqrt(2 - 1250*x*x + 390625*x*x*x*x - 1250*y*y + 390625*y*y*y*y)
         ))/(1248 - 625*Sqrt(2 - 2*x*x + x*x*x*x - 2*y*y + y*y*y*y) + 
       Sqrt(2 - 1250*x*x + 390625*x*x/*exact_solution*/*x*x - 1250*y*y + 390625*y*y*y*y)));
}

double o13(double x, double y)
{
    return 0.35 + 3.25*x*x + 2.5*x*y +
           3.75*y*y +
           Sqrt(0.012500000000000008 +
                7.8125*x*x*x*x + 18.75*x*x*x*y -
                0.5000000000000001*y*y +
                18.75*x*y*y*y + 7.8125*y*y*y*y +
                x*x*(-0.5000000000000001 +
                     46.875*y*y)) -
           1.*Sqrt(Power(-1 + x*x,2)/4. +
                   Power(-0.15000000000000002 +
                         3.75*x*x + 2.5*x*y +
                         3.75*y*y +
                         Sqrt(0.012500000000000008 +
                              7.8125*x*x*x*x +
                              18.75*x*x*x*y -
                              0.5000000000000001*y*y +
                              18.75*x*y*y*y +
                              7.8125*y*y*y*y +
                              x*x*
                              (-0.5000000000000001 +
                               46.875*y*y)),2));
}
double bf13(double x, double y)
{
    return (0.15384615384615383 -
            0.15384615384615383*x*x)/
           (0.10769230769230766 + 1.*x*x +
            0.7692307692307692*x*y +
            1.1538461538461537*y*y +
            0.30769230769230765*
            Sqrt(0.012500000000000008 +
                 7.8125*x*x*x*x +
                 18.75*x*x*x*y -
                 0.5000000000000001*y*y +
                 18.75*x*y*y*y +
                 7.8125*y*y*y*y +
                 x*x*
                 (-0.5000000000000001 +
                  46.875*y*y)));
}

double o12(double x, double y)
{
    return 0.4 + 4.5*x*x + 5.*y*y +
           Sqrt(0.005000000000000003 +
                12.5*x*x*x*x -
                0.5000000000000001*y*y +
                12.5*y*y*y*y +
                x*x*(-0.5000000000000001 +
                     75.*y*y)) -
           1.*Sqrt(Power(-1 + x*x,2)/4. +
                   Power(-0.10000000000000002 +
                         5.*x*x + 5.*y*y +
                         Sqrt(0.005000000000000003 +
                              12.5*x*x*x*x -
                              0.5000000000000001*y*y +
                              12.5*y*y*y*y +
                              x*x*
                              (-0.5000000000000001 +
                               75.*y*y)),2));
}
double o12_2(double x, double y)
{
    return (1 - y*y)/2.;
}
double bf12(double x, double y)
{
    return (0.11111111111111112 -
            0.11111111111111112*x*x)/
           (0.08888888888888888 + 1.*x*x +
            1.1111111111111112*y*y +
            0.22222222222222224*
            Sqrt(0.005000000000000003 +
                 12.5*x*x*x*x -
                 0.5000000000000001*y*y +
                 12.5*y*y*y*y +
                 x*x*
                 (-0.5000000000000001 + 75.*y*y)
                ));
}


double o11(double x, double y)
{
    return  0.25 + x - 1.*y*y -
            1.*Sqrt(x*x -
                    1.*x*(-0.25 + y*y) +
                    1.*Power(-0.25 + y*y,2));
}

double o11_2(double x, double y)
{
    return  2. - x;
}

double bf11(double x, double y)
{
    return (-8.75*(-0.04000000000000001 + y*y)*
            (-0.1 - 0.2857142857142857*x +
             1.*y*y +
             0.2857142857142857*
             Sqrt(x*x -
                  1.*x*(-0.25 + y*y) +
                  1.*Power(-0.25 + y*y,2)) -
             0.2857142857142857*
             Sqrt(6.25*Power(-0.04000000000000001 +
                             y*y,2) -
                  2.5*(-0.04000000000000001 +
                       y*y)*
                  (0.25 + x - 1.*y*y -
                   1.*Sqrt(x*x -
                           1.*x*(-0.25 + y*y) +
                           1.*Power(-0.25 + y*y,2))
                  ) + Power(-0.25 - 1.*x +
                            1.*y*y +
                            Sqrt(x*x -
                                 1.*x*(-0.25 + y*y) +
                                 1.*Power(-0.25 + y*y,2)),
                            2))))/
           (-0.25 - 2.*x + 1.*y*y +
            1.*Sqrt(x*x -
                    1.*x*(-0.25 + y*y) +
                    1.*Power(-0.25 + y*y,2)) -
            1.*Sqrt(x*x -
                    2.5*x*(-0.04000000000000001 +
                           y*y) +
                    6.25*Power(-0.04000000000000001 +
                               y*y,2)) -
            1.*Sqrt(6.25*Power(-0.04000000000000001 +
                               y*y,2) -
                    2.5*(-0.04000000000000001 +
                         y*y)*
                    (0.25 + x - 1.*y*y -
                     1.*Sqrt(x*x -
                             1.*x*(-0.25 + y*y) +
                             1.*Power(-0.25 + y*y,2)))\
                    + Power(-0.25 - 1.*x +
                            1.*y*y +
                            Sqrt(x*x -
                                 1.*x*(-0.25 + y*y) +
                                 1.*Power(-0.25 + y*y,2)),2))
           );
}



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
               (2 - (4*Power(-x + x*x*x,2))/
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
               (2 - (4*Power(-y + y*y*y,2))/
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
                Sqrt(2 - 2*x*x + x*x*x*x/*exact_solution*/ - 2*y*y + y*y*y*y))))/
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
    return 2.*(y*(-8 + 225*x*x + 225*y*y +
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

double f5(double x, double y) {
    return exp(x+y)*2. - 2.*exp(x+y)*((x-X0)*(x-X1)*(y-Y0)*(y-Y1)+(x-X0)*(x-X1)*(y-Y0)+(x-X0)*(x-X1)*(y-Y1)+(x-X0)*(x-X1)+(x-X0)*(y-Y0)*(y-Y1)+(x-X1)*(y-Y0)*(y-Y1)+(y-Y0)*(y-Y1));
}
double u5(double x, double y) {
    return exp(x+y);
}

double f4(double x, double y) {
    return 12.*(y*y*(x*x*x*x-625.) + x*x*(y*y*y*y-1.));
}
double u4(double x, double y) {
    return (x*x*x*x-625.)*(y*y*y*y-1.);
}

double f3(double x, double y) {
    return -2.*sin(x)*sin(y);
}
double u3(double x, double y) {
    return sin(x)*sin(y);
}

double f2(double x, double y) {
    return 2.*(x*x+y*y-2.);
}
double u2(double x, double y) {
    return (x*x-1.)*(y*y-1.);
}

double f1(double x, double y) {
    return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));
}
double u1(double x, double y) {
	//~ printf("done!\n");
    return (x*x*x*x-1.)*(y*y*y*y-1.);
}

double laplace_f(double x, double y) {
    return 0.;
}
double omega_rectangle(double x, double y) {
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}
void init_eq(int id)
{
    omega = 0;
    omega2 = 0;
    structure = 0;
    f_boundary = 0;
    right_part_f = 0;
	u_exact = 0;

    structure = &structure1;
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
        structure = &structure1;
        right_part_f = &f_num;
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
        structure	 = &structure2;
        right_part_f = &f7;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf7;
        omega		 = &o7;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 8)
    {
        structure	 = &structure2;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
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
        structure	 = &structureM;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
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
        structure	 = &structureM;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf10;
        omega		 = &o9;
        omega2		 = &o10_2;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 11)
    {
        structure	 = &structureM;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf11;
        omega		 = &o11;
        omega2		 = &o11_2;
        X0 =   0.;
        X1 =   2.;
        Y0 =  -0.5;
        Y1 =   0.5;
    }
    if(id == 12)
    {
        structure	 = &structureM;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf12;
        omega		 = &o12;
        omega2		 = &o12_2;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 13)
    {   //doesn't work with basises: 1,
        structure	 = &structureM;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf13;
        omega		 = &o13;
        omega2		 = &o12_2;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 14)
    {
        structure = &structure4;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf14;
        omega		 = &o14;
        omega2		 = &o14;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 15)
    {
        structure = &structure1;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf15;
        omega		 = &omega_rectangle;
        omega2		 = &omega_rectangle;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 16)
    {
        structure = &structure4;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf16;
        omega		 = &o16;
        omega2		 = &o16;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    }
    if(id == 17)
    {
        structure = &structure1;
        right_part_f = &f_num;
        //~ u_exact 	 = &u5;
        f_boundary	 = &bf17;
        omega		 = &o17;
        omega2		 = &o17;
        X0 =  -1.;
        X1 =   1.;
        Y0 =  -1.;
        Y1 =   1.;
    } 
}

void GetTask(int id, task *Task)
{
	//~ task returnTask;
    (*Task).Structure = &Structure1;
    if(id == 1)
    {
        (*Task).right_part_f = &f1;
        u_exact 	 	 	 = &u1;
        (*Task).f_boundary	 = &laplace_f;
        (*Task).omega		 = &omega_rectangle;
        (*Task).omega2		 = &omega_rectangle;
        right_part_f = &f1;
        //~ u_exact 	 = &u1;
        f_boundary	 = &laplace_f;
        omega		 = &omega_rectangle;
        omega2		 = &omega_rectangle;
        X0 = -1.;
        X1 =  1.;
        Y0 = -1.;
        Y1 =  1.;        
        (*Task).area.x0 = -1.;
        (*Task).area.x1 =  1.;
        (*Task).area.y0 = -1.;
        (*Task).area.y1 =  1.;
    }

    if(id == 2)
    {
        (*Task).right_part_f  = &f2;
        u_exact 	 	 = &u2;
        (*Task).f_boundary	 = &laplace_f;
        (*Task).omega		 = &omega_rectangle;
        X0 = -1.;
        X1 =  1.;
        Y0 = -1.;
        Y1 =  1.;
        (*Task).area.x0 = -1.;
        (*Task).area.x1 =  1.;
        (*Task).area.y0 = -1.;
        (*Task).area.y1 =  1.;
    }
    //~ return returnTask;
}
