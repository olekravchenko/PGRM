double error_at(double x, double y, task Task)
{
	return fabs(reconstruct_at(Task.solution,x,y)-Task.exact_solution(x,y));
}

double error_C(task Task)
{
	double X, Y, stepx = (Task.area.x1-Task.area.x0)/64., 
				 stepy = (Task.area.y1-Task.area.y0)/64., ret_n = 0.;
	for (X = Task.area.x0; X < Task.area.x1; X+= stepx)
		for (Y = Task.area.y0; Y < Task.area.y1; Y+=stepy)
			ret_n = fmax(ret_n, error_at(X,Y, Task));
		
	return ret_n;
}

double errl1(	basis_args args,
                task *Task,
                int dimension)
{
    double x0 = Task->area.x0;
    double x1 = Task->area.x1;
    double y0 = Task->area.y0;
    double y1 = Task->area.y1;
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
                res += weights[j]*errl1(temp_args, Task, 1);////TODO: complete it
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
                res += weights[j]*error_at(temp_args.x,temp_args.y,*Task);
            }
        }

        return 0.5*res*step;
    }
    return 0.;
}


double errl2(	basis_args args,
				task *Task,
                int dimension)
{
    double x0 = Task->area.x0;
    double x1 = Task->area.x1;
    double y0 = Task->area.y0;
    double y1 = Task->area.y1;
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
                res += weights[j]*errl2(temp_args, Task, 1);////TODO: complete it
            }
        }

        return sqrt(0.5*res*step);
    }
    if (dimension == 1) //case of one dimensional integration
    {
		step = (x1-x0)/intStep;
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < NodesQ; j++)
            {
                temp_args.x = (double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step;
                res += weights[j]*
					error_at(temp_args.x,temp_args.y,*Task)*
					error_at(temp_args.x,temp_args.y,*Task);
            }
        }

        return 0.5*res*step;
    }
    return 0.;
}

