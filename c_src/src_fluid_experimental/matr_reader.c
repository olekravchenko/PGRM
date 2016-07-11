#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
	int n = atoi(argv[1]);
	int i,j;
	double temp;
	FILE *inpt;
	inpt = fopen("matrix","r");
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			fscanf(inpt,"%f",&temp);
			printf("%f\t",temp);
		}
		printf("\n");
	}
	return 0;
}
