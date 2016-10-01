#include <stdio.h>
#include <gsl/gsl_linalg.h>


#define N 10
#define X0 -1.
#define X1  1.
#define Y0 -1.
#define Y1  1.

//depends on which kind of variable supports patricular GPU

//16-bit (tested on CUDA 5.5 + GT8600M + OS X 10.9)
#define diff_step 0.0078125
#define glob_delta 128.

//32-bit (tested on CUDA 5.5 + GTX660 + xUbuntu)
//#define diff_step 0.0009765625
//#define glob_delta 1024.

//64-bit (untested on GPU, tested only on CPU)
//#define diff_step 0.00006103515625
//#define glob_delta 16384.
#define intStep 1.

typedef struct basis_args {
    double	x, y;
    int		m, n;
} basis_args;

typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;

__device__ double phi (double x, double y, int n)
{
    return pow(x,n%N)*pow(y,n/N);
}
__device__ double omega(double x, double y) {
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}
__device__ double structure(double x, double y, int n)
{
    return phi(x,y,n)*omega(x,y);
}

__host__ double Hphi (double x, double y, int n)
{
    return pow(x,n%N)*pow(y,n/N);
}
__host__ double Homega(double x, double y) {
    return (x-X0)*(x-X1)*(y-Y0)*(y-Y1);
}
__host__ double Hstructure(double x, double y, int n)
{
    return Hphi(x,y,n)*Homega(x,y);
}
__host__ double reconstruct_at(gsl_vector *solution,
                               double x, double y)
// Reconstucts value of solution at point (x,y)
{
    int i;
    double result = 0.;
    for(i=0; i<N*N; i++)
        result += gsl_vector_get(solution, i)*Hstructure(x,y,i);

    return result;
}

void plot_region(gsl_vector *solution/*, rect_area plot_area*/)
{
    double hx = (X1-X0)/64.,
           hy = (Y1-Y0)/64.,
           i,j;

    FILE * op;
    op = fopen("../plot_data/plot_region.txt", "w");
    for(i=X0; i<=X1; i+=hx)
        for(j=Y0; j<=Y1; j+=hy)
        {
            fprintf(op, "%15.15f,%15.15f,%15.15f;\n", i,j, reconstruct_at(solution,i,j));
        }
    fclose(op);
    i = system("../bin/plotter.py ../plot_data/plot_region Numerical &");
}


__device__ double right_part_f(double x, double y)
{
    return 12.*(y*y*(x*x*x*x-1.) + x*x*(y*y*y*y-1.));
}


__device__ double left_under_int_new(basis_args arguments)
{
    double 	x = arguments.x;
    double 	y = arguments.y;
    int 	m = arguments.m;
    int 	n = arguments.n;

    return  	structure(x,y,m)*(
                    structure(x+diff_step,y,n)+structure(x-diff_step,y,n)+
                    structure(x,y+diff_step,n)+structure(x,y-diff_step,n)
                    -4.*structure(x,y,n))*glob_delta*glob_delta;
}
__device__ double right_under_int_new(basis_args arguments)
{
    double 	x = arguments.x;
    double 	y = arguments.y;
    int 	m = arguments.m;

    return right_part_f(x,y)*structure(x,y,m);
}


__device__ double nodes[16], weights[16];
__global__ void initGaussInt()
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

__device__ double gauss_integral_right1(rect_area int_area,
                                       basis_args args
                                       )
{
    double x0 = int_area.x0;
    double x1 = int_area.x1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep;

    basis_args temp_args = args;

        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
                res += weights[j]*right_under_int_new(temp_args);
            }
        }

        return 0.5*res*stepx;
}

__device__ double gauss_integral_right2(rect_area int_area,
                                       basis_args args
                                       )
{
    double y0 = int_area.y0;
    double y1 = int_area.y1;

    int i,j;
    double res = 0., stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                res += weights[j]*gauss_integral_right1( int_area, temp_args);
            }
        }

        return 0.5*res*stepy;
}

__device__ double gauss_integral_left1(rect_area int_area,
                                       basis_args args
                                       )
{
    double x0 = int_area.x0;
    double x1 = int_area.x1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep;

    basis_args temp_args = args;

    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 16; j++)
        {
            temp_args.x = (double)(i-1)*stepx + x0 + 0.5*(nodes[j]+1.)*stepx;
            res += weights[j]*left_under_int_new(temp_args);
        }
    }

    return 0.5*res*stepx;
}
__device__ double gauss_integral_left2(rect_area int_area,
                                       basis_args args
                                       )
{
    double y0 = int_area.y0;
    double y1 = int_area.y1;

    int i,j;
    double res = 0., stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 16; j++)
        {
            temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
            res += weights[j]*gauss_integral_left1( int_area, temp_args);
        }
    }

    return 0.5*res*stepy;
}





__global__ void form_matrix_new (float *sys,
                                 float *RightPart)
{
    int i = blockIdx.x, j = threadIdx.x;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
	rect_area int_area = {.x0 = X0, .x1 = X1, .y0 = Y0, .y1 = Y1};
	
    args.m = i;
    
    RightPart[i] = gauss_integral_right2(int_area,args);
    
    args.n = j;
    
    
    sys[i*N*N+j] = gauss_integral_left2(int_area, args);
}


int main()
{
    initGaussInt<<<1,1>>>();

    //as usual, we define pointer to arrays in RAM and GPU RAM
    float *System, *right_part;//, *solution;
    float *dev_System, *dev_right_part;//, *dev_solution;

    //and allocaing this memory
    System = (float *)malloc(N*N * N*N*sizeof(float));
    cudaMalloc (&dev_System, N*N * N*N*sizeof(float));
    right_part = (float *)malloc(N*N*sizeof(float));
    cudaMalloc (&dev_right_part, N*N*sizeof(float));

    //forming the system. TODO: rename this func...
    form_matrix_new<<<N*N, N*N>>>(dev_System, dev_right_part);
    cudaMemcpy( System, dev_System, N*N * N*N*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy( right_part, dev_right_part, N*N*sizeof(float), cudaMemcpyDeviceToHost);


    gsl_matrix 	*Gsys;
    gsl_vector  *Grightpart, *Gsolution;
    Gsys		= gsl_matrix_alloc (N*N,N*N);
    Grightpart	= gsl_vector_alloc(N*N);
    Gsolution	= gsl_vector_alloc(N*N);
    int i, j;

    //filling library-specified mathematical objects
    for(i = 0; i < N*N; i++)
    {
        gsl_vector_set(Grightpart, i, right_part[i]);
        for(j = 0; j < N*N; j++)
        {
            gsl_matrix_set(Gsys, i,j, System[i*N*N+j]);
        }
    }

    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (Gsys, p, &i);
    gsl_linalg_LU_solve (Gsys, p, Grightpart, Gsolution);

	plot_region(Gsolution);


    return 0;
}
