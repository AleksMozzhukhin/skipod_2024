#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h> // Для измерения времени
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N   (2*2*2*2*2*2+2)
float   maxeps = 0.1e-7;
int itmax = 100;
int i,j,k;

float eps;
float A [N][N][N];

void relax();
void init();
void verify();

int main(int an, char **as)
{
    int it;
    clock_t start_time, end_time;
    double total_time, init_time, relax_time, verify_time;

    // Измерение времени функции init
    start_time = clock();
    init();
    end_time = clock();
    init_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Время инициализации: %f секунд\n", init_time);

    // Измерение времени основной итерационной части
    start_time = clock();
    for(it=1; it<=itmax; it++)
    {
        eps = 0.;

        relax();

        if (eps < maxeps) break;
    }
    end_time = clock();
    relax_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Время выполнения relax() в итерациях: %f секунд\n", relax_time);

    // Измерение времени функции verify
    start_time = clock();
    verify();
    end_time = clock();
    verify_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Время выполнения verify(): %f секунд\n", verify_time);

    printf("TOTAL_TIME: %f\n", init_time + relax_time + verify_time);

    return 0;
}


void init()
{
    for(k=0; k<=N-1; k++)
        for(j=0; j<=N-1; j++)
            for(i=0; i<=N-1; i++)
            {
                if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)
                    A[i][j][k]= 0.;
                else
                    A[i][j][k]= ( 4. + i + j + k) ;
            }
}

void relax()
{
    // Первый цикл
    for(k=1; k<=N-2; k++)
        for(j=1; j<=N-2; j++)
            for(i=2; i<=N-3; i++)
            {
                A[i][j][k] = (A[i-1][j][k]+A[i+1][j][k]+A[i-2][j][k]+A[i+2][j][k])/4.;
            }

    // Второй цикл
    for(k=1; k<=N-2; k++)
        for(j=2; j<=N-3; j++)
            for(i=1; i<=N-2; i++)
            {
                A[i][j][k] =(A[i][j-1][k]+A[i][j+1][k]+A[i][j-2][k]+A[i][j+2][k])/4.;
            }

    // Третий цикл
    eps = 0.0; // Инициализируем eps перед циклом
    for(k=2; k<=N-3; k++)
        for(j=1; j<=N-2; j++)
            for(i=1; i<=N-2; i++)
            {
                float e;
                e = A[i][j][k];
                A[i][j][k] = (A[i][j][k-1]+A[i][j][k+1]+A[i][j][k-2]+A[i][j][k+2])/4.;
                eps = Max(eps, fabs(e - A[i][j][k]));
            }
}

void verify()
{
    float s;
    s=0.;

    for(k=0; k<=N-1; k++)
        for(j=0; j<=N-1; j++)
            for(i=0; i<=N-1; i++)
            {
                s = s + A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
            }
    printf("  S = %f\n",s);
}
