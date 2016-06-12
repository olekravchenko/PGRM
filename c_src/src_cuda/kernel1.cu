#include <stdio.h>
//#include <CUDA.h>
//#include <math.h>



#define SIZE	512

__global__ void iterate(float* U, float* Unew)
{
	int m = blockIdx.x*SIZE+threadIdx.x+SIZE+1;
	Unew[m]=0.25*(U[m-SIZE]+U[m+SIZE]+U[m-1]+U[m+1]);
}


int main()
{
	float *U;
	float *dU, *dUnew;
	double h=M_PI/SIZE;
	int n=0,i=0;

	U = (float *)malloc(SIZE*SIZE*sizeof(float));

	cudaMalloc( &dU, 	SIZE*SIZE*sizeof(float));
	cudaMalloc( &dUnew, 	SIZE*SIZE*sizeof(float));
	
	
	
	for(i=0;i<SIZE;i++)
	{
		U[i]=sin(i*h);
		U[(SIZE-1)*SIZE+i]=U[i];
		U[i*SIZE]=-U[i];
		U[i*SIZE+SIZE-1]=-U[i];
	}


	cudaMemcpy( dU, 	U, SIZE*SIZE*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dUnew, 	U, SIZE*SIZE*sizeof(float), cudaMemcpyHostToDevice );



	for(n=0;n<100000;n++){
		iterate<<<SIZE-2,SIZE-2>>>(dU, dUnew);
		cudaMemcpy( dU, dUnew, SIZE*SIZE*sizeof(float), cudaMemcpyDeviceToDevice) ;
	}
	
	
	
	cudaMemcpy( U, dU, SIZE*SIZE*sizeof(float), cudaMemcpyDeviceToHost );


	for(i=0;i<SIZE*SIZE;i+=30)
	{
		printf("%f %f %f\n",h*(i%SIZE),h*(i/SIZE),U[i]);
	}
	
	
	
	free(U);
	cudaFree(dU);
	cudaFree(dUnew);

	return 0;
}
