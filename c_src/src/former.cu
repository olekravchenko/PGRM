#include <stdio.h>

#define N 	10
#define X0 -1.
#define X1  1.
#define Y0 -1.
#define Y1  1.
#define diff_step 0.00000001
#define glob_delta 100000000.
#define intStep 2.

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
    //~ double y0 = int_area.y0;
    //~ double y1 = int_area.y1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep;//, stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
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
    //~ double x0 = int_area.x0;
    //~ double x1 = int_area.x1;
    double y0 = int_area.y0;
    double y1 = int_area.y1;

    int i,j;
    double res = 0.,/* stepx = (x1-x0)/intStep, */stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
        for (i = 1; i <= intStep; i++)
        {
            for (j = 0; j < 16; j++)
            {
                temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
                //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
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
    //~ double y0 = int_area.y0;
    //~ double y1 = int_area.y1;

    int i,j;
    double res = 0., stepx = (x1-x0)/intStep;//, stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
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
    //~ double x0 = int_area.x0;
    //~ double x1 = int_area.x1;
    double y0 = int_area.y0;
    double y1 = int_area.y1;

    int i,j;
    double res = 0.,/* stepx = (x1-x0)/intStep, */stepy = (y1-y0)/intStep;

    basis_args temp_args = args;

    //integral calculations
    for (i = 1; i <= intStep; i++)
    {
        for (j = 0; j < 16; j++)
        {
            temp_args.y = (double)(i-1)*stepy + y0 + 0.5*(nodes[j]+1.)*stepy;
            //res += weights[j]*SubIntegralLeft((*f),x0,x1,(double)(i-1)*step + x0 + 0.5*(nodes[j]+1.)*step,k1,k2);
            res += weights[j]*gauss_integral_left1( int_area, temp_args);
        }
    }

    return 0.5*res*stepy;
}





__global__ void form_matrix_new (float *sys)/*,
                                 //float *RightPart,
                                 rect_area int_area)*/
{
    int i = blockIdx.x, j = threadIdx.x;
    basis_args args;
    args.x = 0.;
    args.y = 0.;
    args.m = 0;
    args.n = 0;
	rect_area int_area = {.x0 = X0, .x1 = X1, .y0 = Y0, .y1 = Y1};
    //~ for(i = 0; i < N*N; i++) // replace this two loops with links to blockId and threadId
    //~ {
    args.m = i;
    
    //RightPart[i] = gauss_integral_right2(int_area,args);
    
    
    //~ for(j = 0; j < N*N; j++)
    //~ {
    args.n = j;
    sys[i*N*N+j] = gauss_integral_left2(int_area, args);
    //~ }
    //~ }
}
//used as example from previous project
__global__ void iter		(float *U, float *Unew, int size)
{
    int k = (blockIdx.x + 1)*size + (threadIdx.x +1);
    float h=0.01;
    Unew[k] = 0.25*(U[k+size]+U[k-size]+U[k-1]+U[k+1]-h*h*2*expf(h*(blockIdx.x+threadIdx.x+2)));
}


int main()
{
    //pointers to host arrays
    float *System;//, *right_part, *solution;
	//rect_area *Area;
    //pointers to device copies of host arrays
    float *dev_System;//, *dev_right_part, *dev_solution;
	
    System = (float *)malloc(N*N * N*N*sizeof(float));
    cudaMalloc( &dev_System, N*N * N*N*sizeof(float));
	int i;//,j;
	for(i = 0; i< N*N * N*N; i++)
	{
		System[i] = 0.;
	}
	cudaMemcpy( dev_System, System, N*N * N*N*sizeof(float), cudaMemcpyHostToDevice);
	form_matrix_new<<<N*N, N*N>>>(dev_System);
	cudaMemcpy( System, dev_System, N*N * N*N*sizeof(float), cudaMemcpyDeviceToHost);
	for(i = 0; i < N*N * N*N; i++)
	{
		printf("%3.3f ",System[i]);
		if(i%(N*N) == 0)
			printf("\n");
	}
	printf("\n");
    return 0;
}
