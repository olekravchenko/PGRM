#include <stdio.h>

int main()
{
	FILE *matr;
	matr =fopen("matrix","r");
	int i,j;
	float cache;
	for(i = 0; i<5; i++)
	{
		for(j = 0; j<5; j++)
		{
			fscanf(matr,"%f", &cache);
			printf("%f\t", cache);
		}
		printf("\n");
	}

		
	return 0;
}