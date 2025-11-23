#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

/* ------------------- STRUCTURES -------------------- */

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

/* PARTITION (venant de ta partie 2 : Tarjan / composantes) */
typedef struct {
    int n;          // nombre de sommets du graphe
    int nb_parts;   // nombre de composantes
    int *part_of;   // part_of[i] = composante du sommet i
} t_partition;

/* ------------------- CREATION -------------------- */

t_matrix createEmptyMatrix(int n);
t_matrix createZeroMatrix(int n);
void freeMatrix(t_matrix *m);

/* ------------------- OPERATIONS -------------------- */

void copyMatrix(const t_matrix *src, t_matrix *dst);
void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C);
int areMultipliable(const t_matrix *A, const t_matrix *B);
double diffMatrices(const t_matrix *M, const t_matrix *N);

/* ------------------- PUISSANCE -------------------- */

t_matrix matrixPower(const t_matrix *M, int k);
void printMatrix(const t_matrix *M, const char *label);

/* ------------------- ADJ -> MATRICE -------------------- */

t_matrix adjacencyListToMatrix(const t_adj_list *list, int n);

/* ------------------- PARTIE 3 -------------------- */

t_matrix subMatrix(t_matrix matrix, t_partition part, int compo_index);
int gcd(int *vals, int nbvals);
int getPeriod(t_matrix sub_matrix);

#endif
