#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TESTS_NUM 4
#define EPS 10e-14
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double antider_left(double x, double h, double m);
double antider_right(double x, double h, double m);
double u(double x, double y);
double norm(double *v1, double *v2, int size);
double **create_EwR(double **A, double w, int size);
void calc_phi(double **R, double **A, double *y, double *b, double *res, double tau, int size);
void upper_triangle_solve(double **A, double *b, double *x, int size);
void lower_triangle_solve(double **A, double *b, double *x, int size);
double **create_A(int size);
double *create_b(double h, int size);
double *alternate_triangular_method(double **A, double *b, double h, int size);
double calc_error(double *y, double h, int size);

double antider_left(double x, double h, double m)
{
    return sin(M_PI * x) / (M_PI * h) - (x / h + 1 - m) * cos(M_PI * x);
}

double antider_right(double x, double h, double m)
{
    return -sin(M_PI * x) / (M_PI * h) - (-x / h + 1 + m) * cos(M_PI * x);
}

double u(double x, double y)
{
    return sin(M_PI * x) * sin(M_PI * y);
}

double norm(double *v1, double *v2, int size)
{
    double res = 0;
    for (int i = 0; i < size; i++) {
        res += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqrt(res);
}

double **create_EwR(double **A, double w, int size)
{
    double **res;
    
    res = (double **)malloc(size * sizeof(double *));
    for (int i = 0; i < size; i++) {
        res[i] = (double *)malloc(size * sizeof(double));
        for (int j = 0; j < size; j++) {
            res[i][j] = w * A[i][j];
        }
        res[i][i] += 1;
    }
    return res;
}

void calc_phi(double **EWR, double **A, double *y, double *b, double *res, double tau, int size) // phi = By - tauAy +taub
{
    double *a = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < size; i++) {
        a[i] = 0;
        for (int j = 0; j <= i; j++) {
            a[i] += EWR[i][j] * y[j];
        }
    }
    for (int i = 0; i < size; i++) {
        res[i] = 0;
        for (int j = i; j < size; j++) {
            res[i] += EWR[i][j] * a[j];
        }
    }
    for (int i = 0; i < size; i++) {
        a[i] = 0;
        for (int j = 0; j < size; j++) {
            a[i] += A[i][j] * y[j];
        }
    }

    for (int i = 0; i < size; i++) {
        res[i] += tau * (b[i] - a[i]); 
    }

    free(a);
}

void upper_triangle_solve(double **A, double *b, double *x, int size)
{
    double s;

    for (int i = size - 1; i > -1; i--) {
        s = 0;
        for (int j = i + 1; j < size; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }
}

void lower_triangle_solve(double **A, double *b, double *x, int size)
{
    double s;

    for (int i = 0; i < size; i++) {
        s = 0;
        for (int j = 0; j < i; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }
}

double **create_A(int size)
{
    double N = size * size;
    double **res = (double **)malloc(N * sizeof(double *));
    int i1, i2, j1, j2;

    for (int i = 0; i < N; i++) {
        res[i] = (double *)malloc(N * sizeof(double));
        i1 = i / size; 
        i2 = i % size;
        for (int j = 0; j < N; j++) {
            j1 = j / size;
            j2 = j % size;
            if (abs(i1 - j1) + abs(i2 - j2) == 0) {
                res[i][j] = 8./3;
            } else if (abs(i1 - j1) + abs(i2 - j2) == 1){
                res[i][j] = 2./3;
            }
        }
    }
    return res;
}

double *create_b(double h, int size)
{
    int i1, i2;
    double *res = (double *)malloc(size * size * sizeof(double));

    for (int i = 0; i < size * size; i++) {
        i1 = i / size + 1;
        i2 = i % size + 1;
        res[i] = (antider_left(i1 * h, h, i1) - antider_left((i1 - 1) * h, h, i1)) + 
                (antider_right((i1 + 1) * h, h, i1) - antider_right(i1 * h, h, i1));
        res[i] *= (antider_left(i2 * h, h, i2) - antider_left((i2 - 1) * h, h, i2)) + 
                (antider_right((i2 + 1) * h, h, i2) - antider_right(i2 * h, h, i2));
        res[i] *= 16 / (3 * M_PI * M_PI * h * h);
    }
    return res;
}

double *alternate_triangular_method(double **A, double *b, double h, int size)
{
    double omega = h / (2 * M_PI);
    double tau =  2 * h / M_PI;
    double **EwR;
    double *phi;
    double *sol_prev, *sol_mid, *sol_next;
    int N = size * size;

    phi = (double *)malloc(N * sizeof(double));
    sol_prev = (double *)malloc(N * sizeof(double));
    sol_mid = (double *)malloc(N * sizeof(double));
    sol_next = (double *)malloc(N * sizeof(double));
    EwR = create_EwR(A, omega, N);

    do {
        for (int i = 0; i < N; i++) {
            sol_prev[i] = sol_next[i];
        }
        calc_phi(EwR, A, sol_prev, b, phi, tau, N);
        upper_triangle_solve(EwR, phi, sol_mid, N);
        lower_triangle_solve(EwR, sol_mid, sol_next, N);
    } while (norm(sol_prev, sol_next, size) > EPS);

    free(phi);
    free(sol_prev);
    free(sol_mid);
    for (int i = 0; i < size; i++) {
        free(EwR[i]);
    }
    free(EwR);
    return sol_next;
}

double calc_error(double *y, double h, int size)
{
    int i1, i2;
    double res = 0;

    for (int i = 0; i < size * size; i++) {
        i1 = i / size;
        i2 = i % size;
        if (fabs(y[i] - u((i1 + 1) * h, (i2 + 1) * h)) > res) {
            res = fabs(y[i] - u((i1 + 1) * h, (i2 + 1) * h));
        }
    }
    return res;
}

int main()
{
    int n;
    double h;
    double **A;
    double *b;
    double *z, *y;
    FILE *fp;

    printf("Enter initial number of points: ");
    scanf("%d", &n);

    h = 1. / (n - 1);
    A = create_A(n - 2);
    b = create_b(h, (n - 2));
    z = alternate_triangular_method(A, b, h, n - 2);

    printf("Error: %.16f\n", calc_error(z, h, n - 2));

    fp = fopen("graph.txt", "w");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ((i == 0) || (i == n - 1) || (j == 0) || (j == n - 1)) {
                fprintf(fp, "%f,%f,%f\n", i * h, j * h, 0.);
            } else {
                fprintf(fp, "%f,%f,%f\n", i * h, j * h, z[(i - 1) * (n - 2) + (j - 1)]);
            }
        }
    }
    fclose(fp);

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(z);

    return 0;
}
