#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define IDX(i, j, N) ((i) * (N) + (j))
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int N;
    int size;
    double h;
    double *x;
    double *b;
    double *r; 
} MGLevel;


double f(double x1, double x2) {
    return 2 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);
}

void allocate_level(MGLevel *level, int N) {
    level->N = N;
    level->size = N * N;
    level->h = 1.0 / (N + 1);
    level->x = (double *)calloc(level->size, sizeof(double));
    level->b = (double *)calloc(level->size, sizeof(double));
    level->r = (double *)calloc(level->size, sizeof(double));
    if (!level->x || !level->b || !level->r) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
}

void free_level(MGLevel *level) {
    free(level->x);
    free(level->b);
    free(level->r);
}

void build_rhs(MGLevel *level) {
    int N = level -> N;
    double h = level -> h;
    for (int i = 0; i < N; ++i) {
        double x1 = (i + 1) * h;
        for (int j = 0; j < N; j++) {
            double x2 = (j + 1) * h;
            level -> b[IDX(i, j, N)]  = f(x1, x2) * h * h;
        }
    }
}

void smooth_jacobi(MGLevel *level, int iterations, double omega) {
    int N = level -> N;
    double *xnew = (double *)calloc(level -> size, sizeof(double));
    for (int it = 0; it < iterations; ++it){
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < N; ++j){
                int idx = IDX(i, j, N);
                double uL = (j > 0) ? level->x[IDX(i, j - 1, N)] : 0;
                double uR = (j < N - 1) ? level->x[IDX(i, j + 1, N)] : 0;
                double uU = (i > 0) ? level->x[IDX(i - 1, j, N)] : 0;
                double uD = (i < N - 1) ? level->x[IDX(i + 1, j, N)] : 0;
                double sigma = uL + uR + uU + uD;
                xnew[idx] = (1 - omega) * level->x[idx] + omega * 0.25 * (level->b[idx] + sigma);
            }
        }
        memcpy(level -> x, xnew, level -> size * sizeof(double));
    }
    free(xnew);
}

void compute_residual(MGLevel *level) {
    int N = level -> N;
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j) {
            int idx = IDX(i, j, N);
            double uC = level -> x[idx];
            double uL = (j > 0) ? level->x[IDX(i, j - 1, N)] : 0;
            double uR = (j < N - 1) ? level->x[IDX(i, j + 1, N)] : 0;
            double uU = (i > 0) ? level->x[IDX(i - 1, j, N)] : 0;
            double uD = (i < N - 1) ? level->x[IDX(i + 1, j, N)] : 0;
            double Ax = 4 * uC - uL - uR - uU - uD;
            level->r[idx] = level->b[idx] - Ax;
        }
    }
}

void restrict_full_weighting(double *fine, double *coarse, int Nf){
    int Nc = Nf / 2;
    for (int i = 0; i< Nc; ++i){
        for(int j = 0; j < Nc; ++j) {
            int fi = 2 * i;
            int fj = 2 * j;
            coarse[IDX(i, j, Nc)] = 0.25 * (fine[IDX(fi, fj, Nf)] +
                                            fine[IDX(fi + 1, fj, Nf)] +
                                            fine[IDX(fi, fj + 1, Nf)] +
                                            fine[IDX(fi + 1, fj + 1, Nf)]);
        }
    }
}

void prolong_linear(double *coarse, double *fine, int Nf){
    int Nc = Nf / 2;
    for (int i = 0; i < Nf; ++i){
        for(int j = 0; j < Nf; ++j) {
            int ci = i/2;
            int cj = j/2;
            fine[IDX(i, j, Nf)] = coarse[IDX(ci, cj, Nc)];
        }
    }
}

void direct_solver(MGLevel *level) {
    smooth_jacobi(level, 50, 0.666);
}

void V_cycle(MGLevel *levels, int l, int lmax, double omega, int nu){
    MGLevel *current = &levels[l];
    smooth_jacobi(current, nu, omega);
    compute_residual(current);
    if(l + 1 <= lmax){
        MGLevel *coarse = &levels[l + 1];
        memset(coarse -> x, 0, coarse -> size * sizeof(double));
        restrict_full_weighting(current -> r, coarse -> b, current -> N);
        if (l + 1 == lmax) {
            direct_solver(coarse);
        } else {
            V_cycle(levels, l + 1, lmax, omega, nu);
        }
        double *e = (double *)calloc(current->size, sizeof(double));
        prolong_linear(coarse -> x, e, current -> N);
        for (int i = 0; i < current -> size; ++i)
            current -> x[i] += e[i];
        free(e);
    }
    smooth_jacobi(current, nu, omega);
}

double residual_norm(MGLevel *level) {
    compute_residual(level);
    double sum = 0.0;
    for (int i = 0; i < level -> size; ++i)
        sum += level -> r[i] * level -> r[i];
    return sqrt(sum);
}

void run_convergence_test(int N, int max_lmax, double omega, int nu, int max_iters){
    double tol = 1e-7;
    printf("Running convergence test for N = %d\n", N);
    for (int lmax = 2; lmax < max_lmax; ++lmax){
        printf("\n --- lmax = %d ---", lmax);
        MGLevel levels[lmax +1];
        for (int l = 0; l <= lmax; ++l){
            allocate_level(&levels[l], N >> l);
            build_rhs(&levels[l]);
        }
        double res = 0.0;

        clock_t start = clock();
        for (int iter = 0; iter < max_iters; ++iter){
            V_cycle(levels, 0, lmax, omega, nu);

            if (iter == 0)
                printf("Initial Residual: %.10e", res);
            if (res <tol){
                printf("Final Residual (Converged at Iteration %d) : %0.10e\n", iter + 1, res);
            }
            if (iter == max_iters)
                printf("Final Residual (Max Iterations Reached): %0.10e\n");
        
        }
        clock_t end = clock();
        printf("Total time elapsed: %.6f seconds", (double)(start-end)/CLOCKS_PER_SEC);
        for(int l =0; l <= lmax; ++l)
            free_level(&levels[l]);
    }
}

int main() {
    int N = 128;
    double omega = 0.666;
    int nu = 3;
    int max_iters = 10000;
    int lmax = 6;
    run_convergence_test(N, lmax, omega, nu, max_iters);
}
