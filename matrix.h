#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct {
    int rows;
    int cols;
    double **data;
} t_matrix;

// Allocation / libération
t_matrix createEmptyMatrix(int n);
t_matrix createZeroMatrix(int n);
void freeMatrix(t_matrix *m);

// Opérations de base 
void copyMatrix(const t_matrix *src, t_matrix *dst);
void multiplyMatrices(const t_matrix *A, const t_matrix *B, t_matrix *C);

// Mesure de différence 
double diffMatrices(const t_matrix *M, const t_matrix *N);

// Puissance entière (M^k)
t_matrix matrixPower(const t_matrix *M, int k);
void printMatrix(const t_matrix *M, const char *label);

/* Construction depuis une liste d’adjacence pondérée (probabilités)
   Format:
   - n sommets
   - pour chaque i: degree d_i, puis d_i couples (j, p_ij) avec j in [0..n-1], p_ij in [0,1]
   Remplit une matrice n x n avec p_ij; les absents valent 0. */
typedef struct {
    int to;
    double prob;
} t_adj_edge;

typedef struct {
    int degree;
    t_adj_edge *edges; // taille degree
} t_adj_list;

t_matrix adjacencyListToMatrix(const t_adj_list *list, int n);

//Outil: vérifie si A et B sont compatibles pour multiplication
int areMultipliable(const t_matrix *A, const t_matrix *B);

#endif
