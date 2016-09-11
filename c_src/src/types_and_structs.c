typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;

typedef struct basis_args {
    double	x, y;
    int		m, n;
} basis_args;

typedef struct task {
    double (*right_part_f)	(double, double);
    double (*f_boundary)	(double, double);
    double (*omega)			(double, double);
    double (*omega2)		(double, double);
    double (*structure)		(double, double, int);
    rect_area area;
    gsl_matrix 	*sys;
    gsl_vector  *rightpart, *solution;
} task;

double Structure1(double x, double y, int n, task Task)
{
    return phi(x,y,n)*Task.omega(x,y);
}
//Neumann problem structure
double Structure2(double x, double y, int n, task Task)
{
    return 	phi(x,y,n)-Task.omega(x,y)*
            ((Task.omega(x+diff_step,y)-Task.omega(x-diff_step,y))*(phi(x+diff_step,y,n)-phi(x-diff_step,y,n))
             +(Task.omega(x,y+diff_step)-Task.omega(x,y-diff_step))*(phi(x,y+diff_step,n)-phi(x,y-diff_step,n)))*
            glob_delta*glob_delta*0.25;
}
//Mixed boundary problem
double StructureM(double x, double y, int n,task Task)
{
    return 	Structure1(x,y,n,Task)-Task.omega(x,y)*Task.omega2(x,y)/(Task.omega(x,y)+Task.omega2(x,y))*
            ((Task.omega2(x+diff_step,y)-Task.omega2(x-diff_step,y))*(Structure1(x+diff_step,y,n,Task)-Structure1(x-diff_step,y,n,Task))
             +(Task.omega2(x,y+diff_step)-Task.omega2(x,y-diff_step))*(Structure1(x,y+diff_step,n,Task)-Structure1(x,y-diff_step,n,Task)))*
            glob_delta*glob_delta*0.25;
}

void tasks_constructor(task *Task, rect_area Area)
{
    Task->sys = gsl_matrix_alloc (N*N,N*N);
    Task->rightpart	= gsl_vector_alloc(N*N);
    Task->solution	= gsl_vector_alloc(N*N);

    Task->area = Area;

    Task->right_part_f = right_part_f;
    Task->f_boundary = f_boundary;
    Task->omega = omega;
    Task->omega2 = omega2;
    Task->structure = structure;


}
