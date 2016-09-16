//double nodes[4], weights[4];

//void initGaussInt()
//{
    ////initializing node & weights
    //nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
    //nodes[0] = -nodes[3];
    //nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
    //nodes[1] = -nodes[2];

    //weights[2] = 0.5+sqrt(30)/36;
    //weights[3] = 0.5-sqrt(30)/36;
    //weights[0] = weights[3];
    //weights[1] = weights[2];
//}
//double nodes[8], weights[8];

//void initGaussInt()
//{
    ////initializing node & weights
    //nodes[0] = -0.1834346424956498;
    //nodes[1] =  0.1834346424956498;
    //nodes[2] = -0.5255324099163290;
    //nodes[3] =  0.5255324099163290;
    //nodes[4] = -0.7966664774136267;
    //nodes[5] =  0.7966664774136267;
    //nodes[6] = -0.9602898564975363;
    //nodes[7] =  0.9602898564975363;

    //weights[0] = 0.3626837833783620;
    //weights[1] = 0.3626837833783620;
    //weights[2] = 0.3137066458778873;
    //weights[3] = 0.3137066458778873;
    //weights[4] = 0.2223810344533745;
    //weights[5] = 0.2223810344533745;
    //weights[6] = 0.1012285362903763;
    //weights[7] = 0.1012285362903763;
//}
double nodes[16], weights[16];

void initGaussInt()
{
    //initializing node & weights
    nodes[0] 	= -0.0950125098376374;
    nodes[1] 	=  0.0950125098376374;
    nodes[2] 	= -0.2816035507792589;
    nodes[3] 	=  0.2816035507792589;
    nodes[4] 	= -0.4580167776572274;
    nodes[5] 	=  0.4580167776572274;
    nodes[6] 	= -0.6178762444026438;
    nodes[7] 	=  0.6178762444026438;
    nodes[8] 	= -0.7554044083550030;
    nodes[9] 	=  0.7554044083550030;
    nodes[10] 	= -0.8656312023878318;
    nodes[11] 	=  0.8656312023878318;
    nodes[12] 	= -0.9445750230732326;
    nodes[13]	=  0.9445750230732326;
    nodes[14] 	= -0.9894009349916499;
    nodes[15] 	=  0.9894009349916499;


    weights[0] 	=  0.1894506104550685;
    weights[1] 	=  0.1894506104550685;
    weights[2] 	=  0.1826034150449236;
    weights[3] 	=  0.1826034150449236;
    weights[4] 	=  0.1691565193950025;
    weights[5] 	=  0.1691565193950025;
    weights[6] 	=  0.1495959888165767;
    weights[7] 	=  0.1495959888165767;
    weights[8] 	=  0.1246289712555339;
    weights[9] 	=  0.1246289712555339;
    weights[10] =  0.0951585116824928;
    weights[11] =  0.0951585116824928;
    weights[12] =  0.0622535239386479;
    weights[13] =  0.0622535239386479;
    weights[14] =  0.0271524594117541;
    weights[15] =  0.0271524594117541;
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

    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
                res += weights[j]*gauss_integral((*f), int_area, temp_args, 1);
            }
        }

        return 0.5*res*stepy;
    }
    if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*(*f)(temp_args);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}
double gauss_integral2(	double (*f)(basis_args,task),
                        rect_area int_area,
                        basis_args args,
                        int dimension,
                        task *Task)
{
    double x0 = int_area.x0;
    double x1 = int_area.x1;
    double y0 = int_area.y0;
    double y1 = int_area.y1;
    //printf("%f %f %f %f\n", x0, x1, y0, y1);
    int i,j;
    double res = 0., stepx = (x1-x0)/intStep, stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
    if (dimension == 2)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
                res += weights[j]*gauss_integral2((*f), int_area, temp_args, 1, Task);
            }
        }

        return 0.5*res*stepy;
    }
    if (dimension == 1)
    {
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*(*f)(temp_args, *Task);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}

