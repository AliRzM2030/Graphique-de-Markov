#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "liste_adjacence.h"

/* ----------------------------------------------------
 *  Conversion List_adj → tableau t_adj_list (pour matrices)
 * ---------------------------------------------------- */

t_adj_list* convertListAdj_to_adjList(List_adj G) {
    t_adj_list *list = malloc(G.taille * sizeof(t_adj_list));

    for (int i = 0; i < G.taille; i++) {
        int degree = 0;
        Cell *c = G.tab[i].head;

        while (c) {
            degree++;
            c = c->next;
        }

        list[i].degree = degree;
        list[i].edges = malloc(degree * sizeof(t_adj_edge));

        c = G.tab[i].head;
        int k = 0;
        while (c) {
            list[i].edges[k].to = c->sommet_d_arrive - 1;
            list[i].edges[k].prob = c->probabilite;
            k++;
            c = c->next;
        }
    }

    return list;
}

/* ---------------------------------------------- */
/*   Allocation interne                            */
/* ---------------------------------------------- */

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

/* ---------------------------------------------- */
/*   Création matrices                             */
/* ---------------------------------------------- */

t_matrix createEmptyMatrix(int n) {
    t_matrix M;
    M.rows = n;
    M.cols = n;
    M.data = alloc2D(n);

    if (!M.data) {
        fprintf(stderr, "Erreur allocation matrice\n");
        M.rows = M.cols = 0;
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
}

/* ---------------------------------------------- */
/* Copie & Multiplication                          */
/* ---------------------------------------------- */

void copyMatrix(const t_matrix *src, t_matrix *dst) {
    if (!src || !dst || !src->data || !dst->data ||
        src->rows != dst->rows || src->cols != dst->cols)
        return;

    memcpy(dst->data[0], src->data[0], src->rows * src->cols * sizeof(double));
}

void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C) {
    memset(C->data[0], 0, C->rows * C->cols * sizeof(double));

    for (int i = 0; i < A->rows; i++) {
        for (int k = 0; k < A->cols; k++) {
            double a = A->data[i][k];
            if (a == 0) continue;

            for (int j = 0; j < B->cols; j++)
                C->data[i][j] += a * B->data[k][j];
        }
    }
}

/* ---------------------------------------------- */
/* Différence                                      */
/* ---------------------------------------------- */

double diffMatrices(const t_matrix *M, const t_matrix *N) {
    double sum = 0;
    int tot = M->rows * M->cols;

    for (int i = 0; i < tot; i++)
        sum += fabs(M->data[0][i] - N->data[0][i]);

    return sum;
}

/* ---------------------------------------------- */
/* Puissance M^k                                   */
/* ---------------------------------------------- */

t_matrix matrixPower(const t_matrix *M, int k) {
    if (!M || !M->data || M->rows != M->cols || k < 0) {
        t_matrix invalid = {0,0,NULL};
        return invalid;
    }

    int n = M->rows;

    if (k == 0) {
        t_matrix I = createZeroMatrix(n);
        for (int i = 0; i < n; i++)
            I.data[i][i] = 1.0;
        return I;
    }

    t_matrix result = createZeroMatrix(n);
    t_matrix base   = createZeroMatrix(n);
    t_matrix tmp    = createZeroMatrix(n);

    for (int i = 0; i < n; i++) result.data[i][i] = 1.0;

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

/* ---------------------------------------------- */
/* Affichage                                       */
/* ---------------------------------------------- */

void printMatrix(const t_matrix *M, const char *label) {
    if (label) printf("%s\n", label);

    for (int i = 0; i < M->rows; i++) {
        for (int j = 0; j < M->cols; j++)
            printf("%.4f ", M->data[i][j]);
        printf("\n");
    }
}

/* ---------------------------------------------- */
/* Conversion adj_list → matrice                   */
/* ---------------------------------------------- */

t_matrix adjacencyListToMatrix(const t_adj_list *list, int n) {
    t_matrix M = createZeroMatrix(n);

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < list[i].degree; k++) {
            int j = list[i].edges[k].to;
            double p = list[i].edges[k].prob;
            M.data[i][j] = p;
        }
    }

    return M;
}

/* ---------------------------------------------- */
/* Sous-matrice + période                          */
/* ---------------------------------------------- */

t_matrix subMatrix(t_matrix matrix, t_partition part, int compo_index) {
    t_classe C = part.classes[compo_index];
    int k = C.nb_vertices;

    t_matrix sub = createZeroMatrix(k);

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            int gi = C.vertices[i] - 1;
            int gj = C.vertices[j] - 1;
            sub.data[i][j] = matrix.data[gi][gj];
        }
    }

    return sub;
}

/* ---------------------------------------------- */
/* GCD + période                                   */
/* ---------------------------------------------- */

int gcd(int *vals, int n) {
    int g = vals[0];

    for (int i = 1; i < n; i++) {
        int a = g, b = vals[i];
        while (b) {
            int t = b;
            b = a % b;
            a = t;
        }
        g = a;
    }
    return g;
}

int getPeriod(t_matrix sub) {
    int n = sub.rows;
    int *vals = malloc(n * sizeof(int));
    int lv = 0;

    t_matrix P = createZeroMatrix(n);
    t_matrix tmp = createZeroMatrix(n);

    copyMatrix(&sub, &P);

    for (int k = 1; k <= n; k++) {
        int diag = 0;
        for (int i = 0; i < n; i++)
            if (P.data[i][i] > 0)
                diag = 1;

        if (diag)
            vals[lv++] = k;

        multiplyMatrices(&P, &sub, &tmp);
        copyMatrix(&tmp, &P);
    }

    int per = (lv > 0 ? gcd(vals, lv) : 0);

    free(vals);
    freeMatrix(&P);
    freeMatrix(&tmp);

    return per;
}