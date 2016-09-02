double nodes[4], weights[4];

void initGaussInt()
{
	//initializing node & weights
	nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
	nodes[0] = -nodes[3];
	nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
	nodes[1] = -nodes[2];
	
	weights[2] = 0.5+sqrt(30)/36;
	weights[3] = 0.5-sqrt(30)/36;
	weights[0] = weights[3];
	weights[1] = weights[2];
}

double gauss_integral(	double (*f)(basis_args),
                        rect_area int_area,
                        basis_args args,
                        int dimension)
{
	double x0 = int_area.x0;
	double x1 = int_area.x1;
	double y0 = int_area.y0;
	double y1 = int_area.y1;
	
    int i,j;
    double res = 0., stepx = (x1-x0)/intStep, stepy = (y1-y0)/intStep;
    
    basis_args temp_args = args;
    //argument processing
	
    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
				temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
                res += weights[j]*gauss_integral((*f), int_area, temp_args, 1);
            }
        }

        return 0.5*res*stepy;
    }
    else if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 4; j++)
            {
				temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*(*f)(temp_args);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}

