#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include "tarjan.h"
#include "liste_adjacence.h"

typedef struct {
    int rows;
    int cols;
    double **data;
} t_matrix;

typedef struct {
    int to;
    double prob;
} t_adj_edge;

typedef struct {
    int degree;
    t_adj_edge *edges;
} t_adj_list;

/* Création */
t_matrix createEmptyMatrix(int n);
t_matrix createZeroMatrix(int n);
void freeMatrix(t_matrix *m);

/* Opérations */
void copyMatrix(const t_matrix *src, t_matrix *dst);
void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C);
int areMultipliable(const t_matrix *A, const t_matrix *B);
double diffMatrices(const t_matrix *M, const t_matrix *N);

/* Puissance */
t_matrix matrixPower(const t_matrix *M, int k);
void printMatrix(const t_matrix *M, const char *label);

/* Conversion List_adj -> Matrix */
t_matrix listAdjToMatrix(List_adj G);

/* PARTIE 3 */
t_matrix subMatrix(t_matrix matrix, t_partition part, int compo_index);
int gcd(int *vals, int nbvals);
int getPeriod(t_matrix sub_matrix);

#endif