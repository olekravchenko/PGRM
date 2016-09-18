

double Structure1(double x, double y, int n, task Task)
{
    return phi(x,y,n)*Task.omega(x,y);
}
double Structure4(double x, double y, int n, task Task)
{
	double temp_omega = Task.omega(x,y);
    return phi(x,y,n)*temp_omega*temp_omega;
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
    if(structure == structure1)
        Task->Structure = Structure1;
    if(structure == structure2)
        Task->Structure = Structure2;
    if(structure == structureM)
        Task->Structure = StructureM;
    if(structure == structure4)
        Task->Structure = Structure4;
	//~ Task->structure = structure;

}
