#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define  Max(a, b) ((a)>(b)?(a):(b))

#define  N   (50)
float maxeps = 0.1e-7;
int itmax = 100;
int i, j, k;

float eps;
float A[N][N][N];

void relax();

void init();

void verify();

int main(int an, char **as) {
    int it;
    double start, end;
    start=omp_get_wtime();
    init();

    for (it = 1; it <= itmax; it++) {
        eps = 0.;

        relax();

        printf("it=%4i   eps=%f\n", it, eps);
        if (eps < maxeps) break;
    }

    verify();
    end=omp_get_wtime();
    printf("TOTAL_TIME: %f\n", end-start);
    return 0;
}


void init() {
#pragma omp parallel for private(i, j, k) collapse(3)
    for (k = 0; k <= N - 1; k++)
        for (j = 0; j <= N - 1; j++)
            for (i = 0; i <= N - 1; i++) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1)
                    A[i][j][k] = 0.;
                else
                    A[i][j][k] = (4. + i + j + k);
            }
}

void relax() {
    eps = 0.0;

#pragma omp parallel
    {
#pragma omp single
        {
            // Первый цикл
            for (k = 1; k <= N - 2; k++)
                for (j = 1; j <= N - 2; j++)
                    for (i = 2; i <= N - 3; i++) {
#pragma omp task firstprivate(i, j, k)
                        {
                            A[i][j][k] = (A[i - 1][j][k] + A[i + 1][j][k] + A[i - 2][j][k] + A[i + 2][j][k]) / 4.;
                        }
                    }
#pragma omp taskwait

            // Второй цикл
            for (k = 1; k <= N - 2; k++)
                for (j = 2; j <= N - 3; j++)
                    for (i = 1; i <= N - 2; i++) {
#pragma omp task firstprivate(i, j, k)
                        {
                            A[i][j][k] = (A[i][j - 1][k] + A[i][j + 1][k] + A[i][j - 2][k] + A[i][j + 2][k]) / 4.;
                        }
                    }
#pragma omp taskwait

            // Третий цикл
            for (k = 2; k <= N - 3; k++)
                for (j = 1; j <= N - 2; j++)
                    for (i = 1; i <= N - 2; i++) {
#pragma omp task firstprivate(i, j, k) shared(eps)
                        {
                            float e;
                            e = A[i][j][k];
                            A[i][j][k] = (A[i][j][k - 1] + A[i][j][k + 1] + A[i][j][k - 2] + A[i][j][k + 2]) / 4.;
                            float local_eps = fabs(e - A[i][j][k]);
#pragma omp critical
                            {
                                eps = Max(eps, local_eps);
                            }
                        }
                    }
#pragma omp taskwait
        }
    }
}

void verify() {
    float s;

    s = 0.;

#pragma omp parallel for private(i, j, k) collapse(3) reduction(+:s)
    for (k = 0; k <= N - 1; k++)
        for (j = 0; j <= N - 1; j++)
            for (i = 0; i <= N - 1; i++) {
                s = s + A[i][j][k] * (i + 1) * (j + 1) * (k + 1) / (N * N * N);
            }
    printf("  S = %f\n", s);
}
