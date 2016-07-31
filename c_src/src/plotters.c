double reconstruct_at(gsl_vector *solution,
                      double x, double y)
// Reconstucts value of solution at point (x,y)
{
    int i;
    double result = 0.;
    for(i=0; i<N*N; i++)
    {
        result += gsl_vector_get(solution, i)*basis(x,y,i);
    }
    return result;
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
            fprintf(op, "%f %f %f\n", i,j, reconstruct_at(solution,i,j));
    fclose(op);

    op = fopen("../plot_data/plot_exact_solution", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%f %f %f\n", i,j, u_exact(i,j));
    fclose(op);

    op = fopen("../plot_data/plot_plot_omega", "w");

    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%f %f %f\n", i,j, omega(i,j));
    fclose(op);

    op = fopen("../plot_data/plot_region_error", "w");
    for(i=x1; i<=x2; i+=hx)
        for(j=y1; j<=y2; j+=hy)
            fprintf(op, "%f %f %f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
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
            fprintf(op, "%f %f %f\n", i,j, reconstruct_at(solution,i,j));
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region Numerical &");
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
            fprintf(op, "%f %f %f\n", i,j, u_exact(i,j));
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
            fprintf(op, "%f %f %f\n", i,j, omega(i,j));
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
            fprintf(op, "%f %f %f\n", i,j, fabs(reconstruct_at(solution,i,j)-u_exact(i,j)));
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region_error Error &");
}
