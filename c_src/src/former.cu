#include <stdio.h>
#include <gsl/gsl_linalg.h>


#define N 10
#define X0 -1.
#define X1  1.
#define Y0 -1.
#define Y1  1.

//depends on which kind of variable supports patricular GPU

//16-bit (tested on CUDA 5.5 + GT8600M + OS X 10.9)
//#define diff_step 0.0078125
//#define glob_delta 128.

//32-bit (tested on CUDA 5.5 + GTX660 + xUbuntu)
//#define diff_step 0.0009765625
//#define glob_delta 1024.

//64-bit (untested on GPU, tested only on CPU)
#define diff_step 0.00006103515625
#define glob_delta 16384.
#define intStep 1.

typedef struct basis_args {
    double	x, y;
    int		m, n;
} basis_args;

typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;


#include "task.cu"



__host__ double reconstruct_at(gsl_vector *solution,
                               double x, double y)
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
    op = fopen("../plot_data/plot_region", "w");
    for(i=X0; i<=X1; i+=hx)
        for(j=Y0; j<=Y1; j+=hy)
            fprintf(op, "%15.15f %15.15f %15.15f\n", i,j, reconstruct_at(solution,i,j));
    fclose(op);
    i = system("screen -d -m ../bin/Plot &");
}


__device__ double left_under_int_new(basis_args arguments)
{
    double 	x = arguments.x;
    double 	y = arguments.y;
    int 	m = arguments.m;
    int 	n = arguments.n;

	//returns S\phi_m \Delta (S\phi_n), where S is operator of the structure
    return			structure(x,y,m)*(
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
//initialization of node & weights values
{
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




/*
 * TODO: try to reduce quantity of similar gauss integral functions, as next to the CPU
 *  way as only possible
 */


__device__ double gauss_integral_right1(rect_area int_area,
                                       basis_args args)
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
                                       basis_args args)
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
                                       basis_args args)
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
                                       basis_args args)
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





__global__ void form_PGRM_sle (float *sys, float *RightPart)
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

__host__ void solve_sle_with_reconstruction(float *System, float *right_part)
{
    gsl_matrix 	*Gsys = gsl_matrix_alloc (N*N,N*N);
    gsl_vector  *Grightpart = gsl_vector_alloc(N*N), 
				*Gsolution = gsl_vector_alloc(N*N);
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
}

__host__ void form_sle_on_gpu (float *System, float *right_part)
{
	float *dev_System, *dev_right_part;

	cudaMalloc (&dev_System, N*N * N*N*sizeof(float));
    cudaMalloc (&dev_right_part, N*N*sizeof(float));

    form_PGRM_sle<<<N*N, N*N>>>(dev_System, dev_right_part);

	cudaMemcpy( System, dev_System, N*N * N*N*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy( right_part, dev_right_part, N*N*sizeof(float), cudaMemcpyDeviceToHost);
	
	cudaFree(dev_right_part);
	cudaFree(dev_System);
}

int main()
{
    initGaussInt<<<1,1>>>();
	//phi = &phi_b3;
    set_b3<<<1,1>>>();
	
	//as usual, we define pointer to arrays in RAM and GPU RAM
    float *System, *right_part;
    System = (float *)malloc(N*N * N*N*sizeof(float));
    right_part = (float *)malloc(N*N*sizeof(float));

	form_sle_on_gpu (System, right_part);
	solve_sle_with_reconstruction(System, right_part);
	
    return 0;
}
