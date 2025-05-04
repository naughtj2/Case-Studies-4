#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IDX(i, j, N) ((i) * (N) + (j))
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double f(double x1, double x2) {
    return 2 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);
}

// Setup for the poisson problem
void poisson_setup(int N, int size, double **A, double *b, double h) {
    *A = (double *) calloc(size * size, sizeof(double));
    if (*A == NULL) {
        perror("Matrix allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        double x1 = (i + 1) * h;
        for (int j = 0; j < N; ++j) {
            double x2 = (j + 1) * h;
            int row = IDX(i, j, N);
            
            b[row] = f(x1, x2);

            (*A)[row * size + row] = 4.0;
            if (i > 0)
                (*A)[row * size + IDX(i - 1, j, N)] = -1.0;
            if (i < N - 1)
                (*A)[row * size + IDX(i + 1, j, N)] = -1.0;
            if (j > 0)
                (*A)[row * size + IDX(i, j - 1, N)] = -1.0;
            if (j < N - 1)
                (*A)[row * size + IDX(i, j + 1, N)] = -1.0;
            
            b[row] *= h * h;
        }

    }
}
int main() {
    int N = 8;
    int size = N * N;
    double h = 1.0 / (N + 1);
    double *A = NULL;
    double *b = (double *)malloc(size * sizeof(double));
    
    if (b == NULL){
        perror("Vector allocation failed");
        return EXIT_FAILURE;
    }

    poisson_setup(N, size, &A, b, h);

    printf("Matrix A:\n");
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j)
            printf("%6.2f ", A[i * size + j]);
        printf("\n");
    }

    printf("\nVector b:\n");
    for (int i = 0; i < size; ++i)
        printf("%6.4f \n", b[i]);
    printf("\n");

    free(A);
    free(b);
    return 0;
}
