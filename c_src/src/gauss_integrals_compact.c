double *nodes, *weights;
int NodesQ;
void initGaussInt(int nodes_quantity)
{
    if(nodes_quantity != 8 && nodes_quantity != 16 && nodes_quantity != 4)
        nodes_quantity = 8;
	NodesQ = nodes_quantity;
        
    init1DArr(&nodes, nodes_quantity);
    init1DArr(&weights, nodes_quantity);


    if(nodes_quantity == 16) {
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
    if(nodes_quantity == 8 ) {
        nodes[0] = -0.1834346424956498;
        nodes[1] =  0.1834346424956498;
        nodes[2] = -0.5255324099163290;
        nodes[3] =  0.5255324099163290;
        nodes[4] = -0.7966664774136267;
        nodes[5] =  0.7966664774136267;
        nodes[6] = -0.9602898564975363;
        nodes[7] =  0.9602898564975363;

        weights[0] = 0.3626837833783620;
        weights[1] = 0.3626837833783620;
        weights[2] = 0.3137066458778873;
        weights[3] = 0.3137066458778873;
        weights[4] = 0.2223810344533745;
        weights[5] = 0.2223810344533745;
        weights[6] = 0.1012285362903763;
        weights[7] = 0.1012285362903763;
    }
    if(nodes_quantity == 4 ) {
        nodes[3] = sqrt(3./7. +2./7.*sqrt(1.2));
        nodes[0] = -nodes[3];
        nodes[2] = sqrt(3./7. -2./7.*sqrt(1.2));
        nodes[1] = -nodes[2];

        weights[2] = 0.5+sqrt(30)/36;
        weights[3] = 0.5-sqrt(30)/36;
        weights[0] = weights[3];
        weights[1] = weights[2];
    }
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
    int i,j;
    double res = 0., step;

    basis_args temp_args = args;

    if (dimension == 2) //case of two dimensional integration in rectangle area
    {
		step = (y1-y0)/intStep;
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < NodesQ; j++)
            {
                temp_args.y = (double)(i-1)*step + y0 + 0.5*(nodes[j]+1.)*step;
                res += /*sqrt(1.-nodes[j]*nodes[j])*/weights[j]*gauss_integral2((*f), int_area, temp_args, 1, Task);
            }
        }

        return 0.5*res*step;
    }
    if (dimension == 1) //case of one dimensional integration
    {
		step = (x1-x0)/intStep;
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < NodesQ; j++)
            {
                temp_args.x = (double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step;
                res += /*sqrt(1.-nodes[j]*nodes[j])*/weights[j]*(*f)(temp_args, *Task);
            }
        }

        return 0.5*res*step;
    }
    return 0.;
}

double gauss_integral3(	double (*f)(basis_args,task),
                        rect_area int_area,
                        basis_args args,
                        int dimension,
                        task *Task)
{
    double x0 = int_area.x0;
    double x1 = int_area.x1;
    double y0 = int_area.y0;
    double y1 = int_area.y1;
    int i,j, M,N;
    double res = 0., stepy = (y1-y0)/intStep, stepx = (x1-x0)/intStep;

    basis_args temp_args = args;

    if (dimension == 2) //case of two dimensional integration in rectangle area
    {
		//step = (y1-y0)/intStep;
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < NodesQ; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
				for (M = 1; M <= intStep; M++)
				{
					for (N = 0; N < NodesQ; N++)
					{
						temp_args.x = (double)(M-1)*stepx + x0 + 0.5*(nodes[N]+1.)*stepx;
						res += weights[j]*weights[N]*(*f)(temp_args, *Task);
					}
				}
            }
        }

        return 0.5*res*stepx*stepy;
    }
    if (dimension == 1) //case of one dimensional integration
    {
		//step = (x1-x0)/intStep;
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < NodesQ; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*(*f)(temp_args, *Task);
            }
        }

        return 0.5*res*stepx;
    }
    return 0.;
}


