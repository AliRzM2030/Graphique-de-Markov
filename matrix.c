#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "liste_adjacence.h"

/* -------------------------- Alloc interne -------------------------- */
static double **alloc2D(int n) {
    if (n <= 0) return NULL;
    double **ptr = malloc(n * sizeof(double *));
    if (!ptr) return NULL;
    ptr[0] = malloc(n * n * sizeof(double));
    if (!ptr[0]) {
        free(ptr);
        return NULL;
    }
    for (int i = 1; i < n; i++)
        ptr[i] = ptr[0] + i * n;
    return ptr;
}

/* -------------------------- Création -------------------------- */
t_matrix createEmptyMatrix(int n) {
    t_matrix M;
    M.rows = n;
    M.cols = n;
    M.data = alloc2D(n);
    if (!M.data) {
        M.rows = M.cols = 0;
        fprintf(stderr, "createEmptyMatrix: alloc failed\n");
    }
    return M;
}

t_matrix createZeroMatrix(int n) {
    t_matrix M = createEmptyMatrix(n);
    if (M.data)
        memset(M.data[0], 0, n * n * sizeof(double));
    return M;
}

void freeMatrix(t_matrix *m) {
    if (!m || !m->data) return;
    free(m->data[0]);
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

/* -------------------------- Utilitaires -------------------------- */
void copyMatrix(const t_matrix *src, t_matrix *dst) {
    if (!src || !dst || !src->data || !dst->data ||
        src->rows != dst->rows || src->cols != dst->cols) {
        fprintf(stderr, "copyMatrix: incompatible matrices\n");
        return;
    }
    memcpy(dst->data[0], src->data[0], src->rows * src->cols * sizeof(double));
}

int areMultipliable(const t_matrix *A, const t_matrix *B) {
    return (A && B && A->data && B->data && A->cols == B->rows);
}

/* -------------------------- Multiplication -------------------------- */
void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C) {
    if (!A || !B || !C || !A->data || !B->data || !C->data) {
        fprintf(stderr, "multiplyMatrices: null pointer\n");
        return;
    }
    if (A->cols != B->rows || C->rows != A->rows || C->cols != B->cols) {
        fprintf(stderr, "multiplyMatrices: dimension mismatch\n");
        return;
    }

    memset(C->data[0], 0, C->rows * C->cols * sizeof(double));
    int n = A->rows;
    int m = A->cols;
    int p = B->cols;

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < m; k++) {
            double a = A->data[i][k];
            if (a == 0.0) continue;
            for (int j = 0; j < p; j++)
                C->data[i][j] += a * B->data[k][j];
        }
    }
}

/* -------------------------- Diff -------------------------- */
double diffMatrices(const t_matrix *M, const t_matrix *N) {
    if (!M || !N || !M->data || !N->data ||
        M->rows != N->rows || M->cols != N->cols) {
        fprintf(stderr, "diffMatrices: incompatible matrices\n");
        return INFINITY;
    }
    double sum = 0;
    int Ntot = M->rows * M->cols;
    for (int i = 0; i < Ntot; i++)
        sum += fabs(M->data[0][i] - N->data[0][i]);
    return sum;
}

/* -------------------------- M^k -------------------------- */
t_matrix matrixPower(const t_matrix *M, int k) {
    t_matrix invalid = {0, 0, NULL};
    if (!M || !M->data || M->rows != M->cols || k < 0)
        return invalid;

    int n = M->rows;
    if (k == 0) {
        t_matrix I = createZeroMatrix(n);
        for (int i = 0; i < n; i++)
            I.data[i][i] = 1.0;
        return I;
    }

    t_matrix result = createZeroMatrix(n);
    t_matrix base = createZeroMatrix(n);
    t_matrix tmp = createZeroMatrix(n);

    for (int i = 0; i < n; i++)
        result.data[i][i] = 1.0; // identite

    copyMatrix(M, &base);

    while (k > 0) {
        if (k & 1) {
            multiplyMatrices(&result, &base, &tmp);
            copyMatrix(&tmp, &result);
        }
        k >>= 1;
        if (k > 0) {
            multiplyMatrices(&base, &base, &tmp);
            copyMatrix(&tmp, &base);
        }
    }

    freeMatrix(&base);
    freeMatrix(&tmp);
    return result;
}

/* -------------------------- Print -------------------------- */
void printMatrix(const t_matrix *M, const char *label) {
    if (label)
        printf("%s\n", label);
    if (!M || !M->data) {
        printf("(null)\n");
        return;
    }
    for (int i = 0; i < M->rows; i++) {
        for (int j = 0; j < M->cols; j++)
            printf("%.4f ", M->data[i][j]);
        printf("\n");
    }
}

/* -------------------------- List_adj → Mat (CORRIGÉ) -------------------------- */
t_matrix listAdjToMatrix(List_adj G) {
    int n = G.taille;
    t_matrix M = createZeroMatrix(n);

    if (!M.data) {
        fprintf(stderr, "listAdjToMatrix: allocation failed\n");
        return M;
    }

    // Parcours de chaque sommet
    for (int i = 0; i < n; i++) {
        Cell *c = G.tab[i].head;
        while (c != NULL) {
            int j = c->sommet_d_arrive - 1;  // CORRECTION: sommets en base 1
            if (j >= 0 && j < n) {
                M.data[i][j] = c->probabilite;
            }
            c = c->next;
        }
    }

    return M;
}

/* -------------------------- Sous-matrice (CORRIGÉ) -------------------------- */
t_matrix subMatrix(t_matrix matrix, t_partition part, int compo_index) {
    if (compo_index < 0 || compo_index >= part.nb_classes) {
        fprintf(stderr, "subMatrix: composante invalide\n");
        return createZeroMatrix(0);
    }

    t_classe C = part.classes[compo_index];
    int k = C.nb_vertices;
    t_matrix sub = createZeroMatrix(k);

    if (!sub.data) return sub;

    // CORRECTION: vertices sont en base 1, il faut convertir en base 0
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            int gi = C.vertices[i] - 1;  // Conversion base 1 → base 0
            int gj = C.vertices[j] - 1;  // Conversion base 1 → base 0

            if (gi >= 0 && gi < matrix.rows && gj >= 0 && gj < matrix.cols) {
                sub.data[i][j] = matrix.data[gi][gj];
            }
        }
    }

    return sub;
}

/* -------------------------- GCD + période -------------------------- */
int gcd(int *vals, int n) {
    if (n == 0) return 0;
    int g = vals[0];
    for (int i = 1; i < n; i++) {
        int a = g, b = vals[i];
        while (b != 0) {
            int tmp = b;
            b = a % b;
            a = tmp;
        }
        g = a;
    }
    return g;
}

int getPeriod(t_matrix sub) {
    int n = sub.rows;
    int *vals = malloc(n * sizeof(int));
    int count = 0;

    t_matrix P = createZeroMatrix(n);
    t_matrix tmp = createZeroMatrix(n);

    copyMatrix(&sub, &P);

    for (int k = 1; k <= n; k++) {
        int diag = 0;
        for (int i = 0; i < n; i++)
            if (P.data[i][i] > 0) {
                diag = 1;
                break;
            }
        if (diag)
            vals[count++] = k;

        multiplyMatrices(&P, &sub, &tmp);
        copyMatrix(&tmp, &P);
    }

    int period = (count > 0) ? gcd(vals, count) : 1;

    free(vals);
    freeMatrix(&P);
    freeMatrix(&tmp);
    return period;
}