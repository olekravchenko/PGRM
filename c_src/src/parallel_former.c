

typedef struct rectangle_area {
    double x1,x2;
    double y1,y2;
} rectangle_area;

typedef struct parallel_arg {
    rectangle_area area;
    int i, j;
    double result;
} parallel_arg;

typedef struct parallel_arg2 {
    rectangle_area area;
    int threadId;
    double **result;
} parallel_arg2;

void *parallel_pre_former(void *arg)
{
    parallel_arg *a = (parallel_arg *)arg;
    a->result = integralLeft(left_under_int,a->area.x1,a->area.x2,a->area.y1,a->area.y2,a->i,a->j);
    return 0;
}

void *parallel_pre_former2(void *arg)
{
    parallel_arg2 *a = (parallel_arg2 *)arg;
    int id = a->threadId;
    int i, j;
    for(i = 0; i<N*N/16; i++)
    {
        for(j = 0; j<N*N; j++)
            a->result[i+id*N*N/16][j] = integralLeft(left_under_int,a->area.x1,a->area.x2,a->area.y1,a->area.y2,i+id*N*N/16,j);
    }
    return 0;
}

/*
 * todo: fix crashing on non 2^n basis function launch
 * 
 */
void form_matrix_parallel(gsl_matrix * system,
                          gsl_vector * RightPart,
                          double x1, double x2,
                          double y1, double y2)
// Forms SLE system
// system 	- left part matrix form of system
// RightPart- right part vector of coefficients
// x1, x2	- sizes of rectangle by x
// y1, y2	- sizes of rectangle by y
{
    int i, j, k, return_code;
    parallel_arg2 arg[16];
    rectangle_area general_area;
    general_area.x1 = x1;
    general_area.x2 = x2;
    general_area.y1 = y1;
    general_area.y2 = y2;

    double **prematrixarray;
    init2DArr(&prematrixarray, N*N, N*N);
    pthread_t threads[16];

    for(i = 0; i < 16; i++)
    {
        arg[i].area = general_area;
        arg[i].threadId = i;
        arg[i].result = prematrixarray;
        return_code=pthread_create(&threads[i],NULL,parallel_pre_former2, (void*)&arg[i]);
    }

    for(k = 0; k < 16; k++)
    {
        return_code = pthread_join(threads[k], NULL);

        for(i = 0; i<N*N/16; i++)
            for(j = 0; j < N*N; j++)
                prematrixarray[i+k*N*N/16][j]=arg[i].result[i+k*N*N/16][j];
    }
    return_code++;
    for(i = 0; i < N*N; i++)
    {
		//for non variational solver remove minus
        gsl_vector_set(RightPart, i, integralRight(right_part_f,basis,x1,x2,y1,y2,i));
        for(j = 0; j < N*N; j++)
            gsl_matrix_set(system, i,j, prematrixarray[i][j]);
    }
    free2DArr(&prematrixarray, N*N);
}
