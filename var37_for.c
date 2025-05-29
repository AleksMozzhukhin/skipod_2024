#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h> // Подключаем OpenMP
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N   (2*2*2*2*2*2+2)
float   maxeps = 0.1e-7;
int itmax = 100;

// Объявляем переменную eps как глобальную, но она будет использоваться в редукции
float eps;
float A [N][N][N];

void relax();
void init();
void verify();

int main(int an, char **as)
{
    int it;
    double start_time, end_time;
    double init_time, relax_time, verify_time;
    double program_start_time, program_end_time;
    double total_time;

    program_start_time = omp_get_wtime(); // Начало общего времени выполнения

    // Измерение времени функции init
    start_time = omp_get_wtime();
    init();
    end_time = omp_get_wtime();
    init_time = end_time - start_time;
    printf("Время инициализации: %f секунд\n", init_time);

    // Измерение времени основной итерационной части
    start_time = omp_get_wtime();
    for(it=1; it<=itmax; it++)
    {
        eps = 0.0; // Инициализируем eps перед каждой итерацией

        relax();

        if (eps < maxeps) break;
    }
    end_time = omp_get_wtime();
    relax_time = end_time - start_time;
    printf("Время выполнения relax() в итерациях: %f секунд\n", relax_time);

    // Измерение времени функции verify
    start_time = omp_get_wtime();
    verify();
    end_time = omp_get_wtime();
    verify_time = end_time - start_time;
    printf("Время выполнения verify(): %f секунд\n", verify_time);

    total_time = init_time + relax_time + verify_time;
    printf("TOTAL_TIME: %f\n", total_time);
    return 0;
}

void init()
{
    // Параллелизуем инициализацию массива
#pragma omp parallel for collapse(3) shared(A)
    for(int k=0; k<=N-1; k++)
        for(int j=0; j<=N-1; j++)
            for(int i=0; i<=N-1; i++)
            {
                if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)
                    A[i][j][k]= 0.;
                else
                    A[i][j][k]= ( 4.0 + i + j + k );
            }
}

void relax()
{
    // Первый цикл: обновление по оси i
#pragma omp parallel for collapse(3) shared(A)
    for(int k=1; k<=N-2; k++)
        for(int j=1; j<=N-2; j++)
            for(int i=2; i<=N-3; i++)
            {
                A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i-2][j][k] + A[i+2][j][k]) / 4.0;
            }

    // Второй цикл: обновление по оси j
#pragma omp parallel for collapse(3) shared(A)
    for(int k=1; k<=N-2; k++)
        for(int j=2; j<=N-3; j++)
            for(int i=1; i<=N-2; i++)
            {
                A[i][j][k] = (A[i][j-1][k] + A[i][j+1][k] + A[i][j-2][k] + A[i][j+2][k]) / 4.0;
            }

    // Третий цикл: обновление по оси k и вычисление eps
    // Используем редукцию для eps
#pragma omp parallel for collapse(3) shared(A) reduction(max:eps)
    for(int k=2; k<=N-3; k++)
        for(int j=1; j<=N-2; j++)
            for(int i=1; i<=N-2; i++)
            {
                float e = A[i][j][k];
                A[i][j][k] = (A[i][j][k-1] + A[i][j][k+1] + A[i][j][k-2] + A[i][j][k+2]) / 4.0;
                float local_eps = fabs(e - A[i][j][k]);
                eps = Max(eps, local_eps); // OpenMP автоматически выполняет редукцию
            }
}

void verify()
{
    float s = 0.0;

    // Параллелизуем вычисление суммы
#pragma omp parallel for collapse(3) reduction(+:s) shared(A)
    for(int k=0; k<=N-1; k++)
        for(int j=0; j<=N-1; j++)
            for(int i=0; i<=N-1; i++)
            {
                s += A[i][j][k] * (i+1) * (j+1) * (k+1) / (float)(N * N * N);
            }
    printf("  S = %f\n", s);
}
