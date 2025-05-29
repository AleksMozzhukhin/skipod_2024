#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#define  Max(a,b) ((a)>(b)?(a):(b))
#define  Min(a,b) ((a)<(b)?(a):(b))
#define  N   (2*2*2*2*2*2+2)
#define TAG_DOWN 0
#define TAG_UP   1

float maxeps = 0.1e-7;
int itmax = 100;

float eps;
float A_local[N][N][N];

float halo_down[2][N][N];
float halo_up[2][N][N];

void relax(int rank, int size, int k_start, int k_end);
void init(int rank, int size, int k_start, int k_end);
void verify(int rank, int size, int k_start, int k_end);
void exchange_halos(int rank, int size, int k_start, int k_end);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int base = N / size;
    int rem = N % size;
    int local_K;
    int k_start, k_end;
    if (rank < rem) {
        local_K = base + 1;
        k_start = rank * local_K;
    } else {
        local_K = base;
        k_start = rem * (base + 1) + (rank - rem) * base;
    }
    k_end = k_start + local_K - 1;

    double program_end_time, program_start_time = omp_get_wtime();

    init(rank, size, k_start, k_end);

    for(int it = 1; it <= itmax; it++)
    {
        eps = 0.0f;
        relax(rank, size, k_start, k_end);

        float global_eps;
        MPI_Allreduce(&eps, &global_eps, 1, MPI_FLOAT,
        MPI_MAX, MPI_COMM_WORLD);
        eps = global_eps;
        if (eps < maxeps) break;
    }

    verify(rank, size, k_start, k_end);

    program_end_time = omp_get_wtime();
    if(rank == 0) printf("TOTAL_TIME: %f\n",
    program_end_time - program_start_time);

    MPI_Finalize();
    return 0;
}

void init(int rank, int size, int k_start, int k_end)
{
    for(int kk = k_start; kk <= k_end; kk++) {
        int k_local = kk - k_start;

        if (k_local < 0 || k_local >= N) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for(int j = 0; j < N; j++) {
            for(int i = 0; i < N; i++) {
                if(i == 0 || i == N-1 || j == 0 ||
                j == N-1 || kk == 0 || kk == N-1)
                    A_local[i][j][k_local] = 0.0f;
                else
                    A_local[i][j][k_local] = 4.0f + i + j + kk;
            }
        }
    }

}

void exchange_halos(int rank, int size, int k_start, int k_end)
{
    MPI_Status status;
    int left = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int right = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    int local_K = k_end - k_start + 1;

    float send_down[2][N][N];
    float send_up[2][N][N];

    for(int j = 0; j < N; j++) {
        for(int i = 0; i < N; i++) {
            if(local_K >=2){
                send_down[0][i][j] = A_local[i][j][0];
                send_down[1][i][j] = A_local[i][j][1];

                send_up[0][i][j] = A_local[i][j][local_K - 2];
                send_up[1][i][j] = A_local[i][j][local_K - 1];
            }
            else{
                send_down[0][i][j] = 0.0f;
                send_down[1][i][j] = 0.0f;
                send_up[0][i][j] = 0.0f;
                send_up[1][i][j] = 0.0f;
            }
        }
    }

    MPI_Sendrecv(send_down, 2 * N * N, MPI_FLOAT, left, TAG_DOWN,
                 halo_up, 2 * N * N, MPI_FLOAT, left, TAG_UP,
                 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(send_up, 2 * N * N, MPI_FLOAT, right, TAG_UP,
                 halo_down, 2 * N * N, MPI_FLOAT, right, TAG_DOWN,
                 MPI_COMM_WORLD, &status);

}

void relax(int rank, int size, int k_start, int k_end)
{
    int local_K = k_end - k_start +1;

    for(int kk = Max(k_start,1); kk <= Min(k_end, N-2); kk++) {
        if(local_K >=1){
            int k_local = kk - k_start;
            for(int j = 1; j <= N-2; j++) {
                for(int i = 2; i <= N-3; i++) {
                    A_local[i][j][k_local] =
                    (A_local[i-1][j][k_local] +
                    A_local[i+1][j][k_local] +
                    A_local[i-2][j][k_local] +
                    A_local[i+2][j][k_local]) /
                    4.0f;
                }
            }
        }
    }
    for(int kk = Max(k_start,1); kk <= Min(k_end, N-2); kk++) {
        if(local_K >=1){
            int k_local = kk - k_start;
            for(int j = 2; j <= N-3; j++) {
                for(int i = 1; i <= N-2; i++) {
                    A_local[i][j][k_local] =
                    (A_local[i][j-1][k_local] +
                    A_local[i][j+1][k_local] +
                    A_local[i][j-2][k_local] +
                    A_local[i][j+2][k_local])
                    / 4.0f;
                }
            }
        }
    }

    if(local_K >=5){
        exchange_halos(rank, size, k_start, k_end);

        for(int kk = Max(k_start,2); kk <= Min(k_end, N-3); kk++) {
            int k_local = kk - k_start;
            if(k_local <2 || k_local > (local_K -3)){
                continue;
            }

            for(int j =1; j <=N-2; j++) {
                for(int i=1; i <=N-2; i++) {
                    float e = A_local[i][j][k_local];
                    float a_km1 = (k_local -1 >=0) ?
                    A_local[i][j][k_local -1] : halo_down[0][i][j];
                    float a_km2 = (k_local -2 >=0) ?
                    A_local[i][j][k_local -2] : halo_down[1][i][j];
                    float a_kp1 = (k_local +1 < local_K) ?
                    A_local[i][j][k_local +1] : halo_up[0][i][j];
                    float a_kp2 = (k_local +2 < local_K) ?
                    A_local[i][j][k_local +2] : halo_up[1][i][j];
                    A_local[i][j][k_local] =
                    (a_km1 + a_kp1 + a_km2 + a_kp2) / 4.0f;
                    float diff = fabsf(e - A_local[i][j][k_local]);
                    eps = Max(eps, diff);
                }
            }
        }
    }
}

void verify(int rank, int size, int k_start, int k_end)
{
    float s_local = 0.0f;
    int local_K = k_end - k_start +1;

    for(int kk = k_start; kk <= k_end; kk++) {
        if(local_K >=1){
            int k_local = kk - k_start;
            if (k_local <0 || k_local >= N){
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            for(int j = 0; j < N; j++) {
                for(int i = 0; i < N; i++) {
                    if(i <0 || i >=N || j <0 || j >=N){
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }

                    s_local += A_local[i][j][k_local] * (i + 1)
                    * (j + 1) * (kk + 1) / (float)(N * N * N);
                }
            }
        }
    }

    float s;
    MPI_Reduce(&s_local, &s, 1, MPI_FLOAT, MPI_SUM,
            0, MPI_COMM_WORLD);

    if(rank == 0) {
        printf("  S = %f\n", s);
    }
}
