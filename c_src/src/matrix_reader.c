#include <stdio.h>

int main()
{
	FILE *matr;
	matr =fopen("matrix","r");
	int i,j;
	float cache;
	for(i = 0; i<9; i++)
	{
		for(j = 0; j<9; j++)
		{
			fscanf(matr,"%f", &cache);
			printf("%3.3g\t", cache);
		}
		printf("\n");
	}

		
	return 0;
}
