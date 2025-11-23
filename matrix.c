#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

static double **alloc2D(int n) {
    double **ptr = (double **)malloc(n * sizeof(double *));
    if (!ptr) return NULL;
    ptr[0] = (double *)malloc(n * n * sizeof(double));
    if (!ptr[0]) { free(ptr); return NULL; }
    for (int i = 1; i < n; i++) {
        ptr[i] = ptr[0] + i * n;
    }
    return ptr;
}

t_matrix createEmptyMatrix(int n) { // cree une matrice n x n avec les probas de passage
    t_matrix M;
    M.rows = n; M.cols = n;
    M.data = alloc2D(n);
    if (!M.data) {
        fprintf(stderr, "Allocation failed in createEmptyMatrix\n");
        M.rows = M.cols = 0;
    }
    return M;
}

t_matrix createMatrixZero (int n) {
    t_matrix M = createEmptyMatrix(n);
    if (M.data) {
        memset(M.data[0], 0, n * n * sizeof(double));
    }
    return M;
}

void freeMatrix(t_matrix *m) {
    if (!m || !m->data) return;
    free(m->data[0]);
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

void copyMatrix(const t_matrix *src, t_matrix *dst) {
    if (!src || !dst || !src->data || !dst->data || src->rows != dst->rows || src->cols != dst->cols) {
        fprintf(stderr, "copyMatrix: incompatible matrices\n");
        return;
    }
    memcpy(dst->data[0], src->data[0], src->rows * src->cols * sizeof(double));
}

int areMultipliable(const t_matrix *A, const t_matrix *B) {
    return A && B && A->cols == B->rows;
}

void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C) {
    if (!A || !B || !C || !A->data || !B->data || !C->data) {
        fprintf(stderr, "multiplyMatrices: null matrix pointer\n");
        return;
    }
    if (A->cols != B->rows || C->rows != A->rows || C->cols != B->cols) {
        fprintf(stderr, "multiplyMatrices: dimension mismatch\n");
        return;
    }

    int n = A->rows;
    int m = A->cols; 
    int p = B->cols;

    // Triple boucle classique
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < m; k++) {
            double a_ik = A->data[i][k];
            if (a_ik == 0.0) continue;
            for (int j = 0; j < p; j++) {
                C->data[i][j] += a_ik * B->data[k][j];
            }
        }
    }
}

double diffMatrices(const t_matrix *M, const t_matrix *N) {
    if (!M || !N || !M->data || !N->data || M->rows != N->rows || M->cols != N->cols) {
        fprintf(stderr, "diffMatrices: incompatible matrices\n");
        return INFINITY;
    }
    double s = 0.0;
    int n = M->rows * M->cols;
    for (int t = 0; t < n; t++) {
        s += fabs(M->data[0][t] - N->data[0][t]);
    }
    return s;
}

t_matrix matrixPower(const t_matrix *M, int k) {
    if (!M || !M->data || M->rows != M->cols || k < 0) {
        fprintf(stderr, "matrixPower: invalid input\n");
        return createZeroMatrix(0);
    }
    int n = M->rows;

    // Identité pour k=0
    if (k == 0) {
        t_matrix I = createZeroMatrix(n);
        for (int i = 0; i < n; i++) I.data[i][i] = 1.0;
        return I;
    }

    // Exponentiation rapide (binary exponentiation)
    t_matrix result = createZeroMatrix(n);
    for (int i = 0; i < n; i++) result.data[i][i] = 1.0;

    t_matrix base = createEmptyMatrix(n);
    copyMatrix(M, &base);

    t_matrix tmp = createEmptyMatrix(n);

    int e = k;
    while (e > 0) {
        if (e & 1) {
            multiplyMatrices(&result, &base, &tmp);
            copyMatrix(&tmp, &result);
        }
        e >>= 1;
        if (e > 0) {
            multiplyMatrices(&base, &base, &tmp);
            copyMatrix(&tmp, &base);
        }
    }

    freeMatrix(&base);
    freeMatrix(&tmp);
    return result;
}

void printMatrix(const t_matrix *M, const char *label) {
    if (label) printf("%s\n", label);
    if (!M || !M->data) { printf("(null)\n"); return; }
    for (int i = 0; i < M->rows; i++) {
        for (int j = 0; j < M->cols; j++) {
            printf("%.2f%s", M->data[i][j], (j + 1 == M->cols) ? "" : " ");
        }
        printf("\n");
    }
}

t_matrix adjacencyListToMatrix(const t_adj_list *list, int n) {
    t_matrix M = createZeroMatrix(n);
    if (!list || !M.data) return M;
    for (int i = 0; i < n; i++) {
        int d = list[i].degree;
        for (int k = 0; k < d; k++) {
            int j = list[i].edges[k].to;
            double p = list[i].edges[k].prob;
            if (j < 0 || j >= n) {
                fprintf(stderr, "adjacencyListToMatrix: invalid edge index %d\n", j);
                continue;
            }
            M.data[i][j] = p;
        }
    }
    return M;
}
