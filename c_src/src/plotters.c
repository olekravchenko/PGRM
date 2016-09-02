double reconstruct_at(gsl_vector *solution,
                      double x, double y)
// Reconstucts value of solution at point (x,y)
{
    int i;
    double result = 0.;

    if(x==X0 && y==Y0)
        return reconstruct_at(solution, x+5.*diff_step,y+5.*diff_step);
    else if(x==X0 && y==Y1)
        return reconstruct_at(solution, x+5.*diff_step,y-5.*diff_step);
    else if(x==X1 && y==Y0)
        return reconstruct_at(solution, x-5.*diff_step,y+5.*diff_step);
    else if(x==X1 && y==Y1)
        return reconstruct_at(solution, x-5.*diff_step,y-5.*diff_step);
    for(i=0; i<N*N; i++)
    {
        result += gsl_vector_get(solution, i)*structure(x,y,i);
    }

    result+=f_boundary(x,y);
    return result;

}

double DeltaU(gsl_vector *solution,double x, double y)
{
	return (reconstruct_at(solution,x+diff_step,y)+reconstruct_at(solution,x-diff_step,y)+
			reconstruct_at(solution,x,y+diff_step)+reconstruct_at(solution,x,y-diff_step)-
			4.*reconstruct_at(solution,x,y))*glob_delta*glob_delta;
}

void multiplot(gsl_vector *solution,
               double x1, double x2,
               double y1, double y2)
{
    double hx = (x2-x1)/64.,
           hy = (y2-y1)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, reconstruct_at(solution,i,j));
    fclose(op);

    op = fopen("../plot_data/plot_exact_solution", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, u_exact(i,j));
    fclose(op);

    op = fopen("../plot_data/plot_plot_omega", "w");

    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, omega(i,j));
    fclose(op);

    op = fopen("../plot_data/plot_region_error", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
    fclose(op);

    i = system("../bin/plotter_experimental.py &");
}


void plot_region(gsl_vector *solution,
                 double x1, double x2,
                 double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
    double hx = (x2-x1)/64.,
           hy = (y2-y1)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
        {
				fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, reconstruct_at(solution,i,j));
		}
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region Numerical &");
}
void plot_region_colorplot(gsl_vector *solution,
                           double x1, double x2,
                           double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
    double hx = (x2-x1)/64.,
           hy = (y2-y1)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
        {
				fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, reconstruct_at(solution,i,j));
		}
    fclose(op);
    i = system("../bin/plotter_colorplot.py ../plot_data/plot_region &");
}
void plot_laplacian(gsl_vector *solution,
                 double x1, double x2,
                 double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
    double hx = (x2-x1)/64.,
           hy = (y2-y1)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_laplacian", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
        {
				fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, DeltaU(solution,i,j));
		}
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_laplacian Numerical &");
}

void plot_exact_solution(double x1, double x2,
                         double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
    double	hx = (x2-x1)/64.,
            hy = (y2-y1)/64.,
            i,j;

    FILE * op;
    op = fopen("../plot_data/plot_exact_solution", "w");

    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, u_exact(i,j));
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_exact_solution Exact &");
}

void plot_omega
(double x1, double x2,
 double y1, double y2)
// Plot solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y

{
    double	hx = (x2-x1)/64.,
            hy = (y2-y1)/64.,
            i,j;

    FILE * op;
    op = fopen("../plot_data/plot_plot_omega", "w");

    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, omega(i,j));
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_plot_omega Omega &");
}

void plot_region_error(gsl_vector *solution,
                       double x1, double x2,
                       double y1, double y2)
// Plot abs error of solution in rectangle region
// from x1 till x2 by x, and from y1 till y2 by y
{
    double hx = (x2-x1)/64.,
           hy = (y2-y1)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region_error", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region_error Error &");
}
