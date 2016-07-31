#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*	polynomial-calculated atomic functions						*/
/*	version 0.1 + poly										*/
/*															*/
/*	structure: 												*/
/*	f_XXXX(double x, ... ) - atomic function "XXXX"			*/
/*	F_XXXX(double t, ... ) - image of "XXXX"				*/
/*	o_XXXX(...) - any other function, used further			*/
/*	f_XXXX_poly - atomic function "XXXX" in polynomial form	*/
/*															*/

double f_up_poly(double x)
{
	x = fabs(x);
	if(0.875 < x && x <= 1.)
		return 64./3 *(x-1.)*(x-1.)*(x-1.)*(x-1.);
	if(0.625 < x && x <= 0.875)
		return 0.010416667*(-353.+2784.*x-6528.*x*x+6144.*x*x*x-2048.*x*x*x*x);
	if(0.5 < x && x <= 0.625)
		return 0.166666667*(17.-76.*x+192.*x*x-256.*x*x*x+128.*x*x*x*x);
	if(0.375 < x && x <= 0.5)
		return 0.166666667*(1.+52.*x-192.*x*x+256.*x*x*x-128.*x*x*x*x);
	if(0.125 < x && x <= 0.375)
		return 0.010416667*(97.-32.*x+384.*x*x-2048.*x*x*x+2048.*x*x*x*x);
	if(0.< x && x <= 0.125)
		return 0.33333333333*(3.-64.*x*x*x*x);
	return 0.;
}

double f_fup3_poly(double x)
{
	x = fabs(x);
	if(2.75 <= x && x < 3.)
		return -(0.006349206)*pow((-3 + x),7);
	if(1.25 <= x && x <= 1.75)
		return 0.8044712611607143 - (10937.*x)/9216. + (543.*pow(x,2))/256. - (2005.*x*x*x)/576. + 
   (49.*pow(x,4))/16. - (17.*pow(x,5))/12. + pow(x,6)/3. - (2.*pow(x,7))/63.;
	if(0 <= x && x <= 0.25)
		return 0.5800130208333333 - (37.*x*x)/64. + pow(x,4)/4. - (4.*pow(x,7))/63.;
	if(0.75 <= x && x <= 1.)
		return 0.5969556051587301 - (91.*x)/576. + (5.*x*x)/96. - (25.*x*x*x)/18. + 
   (37.*pow(x,4))/18. - (4.*pow(x,5))/3. + (4.*pow(x,6))/9. - (4.*pow(x,7))/63.;
	if(2. <= x && x <= 2.25)
		return (1605851. - 3540208.*x + 3765552.*x*x - 2562560.*x*x*x + 1173760.*pow(x,4) - 
     344064.*pow(x,5) + 57344.*pow(x,6) - 4096.*pow(x,7))/645120.;
	if(1. < x && x < 1.25)
		return 0.5017175099206349 + (293.*x)/576. - (187.*x*x)/96. + (35.*x*x*x)/18. - 
   (23.*pow(x,4))/18. + (2.*pow(x,5))/3. - (2.*pow(x,6))/9. + (2.*pow(x,7))/63.;
	if(2.25 < x && x < 2.75)
		return (-1571267. + 7799932.*x - 12309360.*x*x + 9571520.*x*x*x- 4184320.*x*x*x*x + 
 1053696.*pow(x,5)- 143360.*pow(x,6) + 8192.*pow(x,7))/1290240.;
	if(1.75 < x && x < 2.)
		return -2.386962115575397 + (13339.*x)/1152. - (7589.*x*x)/384. + (625.*x*x*x)/36. - 
   (637.*pow(x,4))/72. + (8.*pow(x,5))/3. - (4.*pow(x,6))/9. + (2.*pow(x,7))/63.;
	if(0.25 < x && x < 0.75)
		return 0.5800052703373015 + x/4608. - (223.*x*x)/384. + (5.*x*x*x)/288. + 
   (13.*pow(x,4))/72. + pow(x,5)/6. - (2.*pow(x,6))/9. + (4.*pow(x,7))/63.;
	return 0.;
}

