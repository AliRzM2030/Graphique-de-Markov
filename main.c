#include <stdio.h>
#include "liste_adjacence.h"
#include "tarjan.h"
#include "hasse.h"
#include "matrix.h"

void run_part3_demo(void) {
    const int n = 5;
    t_matrix M = createZeroMatrix(n);

    double vals[5][5] = {
        {0.34, 0.27, 0.00, 0.18, 0.21},
        {0.20, 0.40, 0.20, 0.00, 0.20},
        {0.00, 0.41, 0.37, 0.09, 0.13},
        {0.00, 0.68, 0.20, 0.12, 0.00},
        {0.12, 0.30, 0.00, 0.00, 0.58}
    };

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            M.data[i][j] = vals[i][j];

    printMatrix(&M, "\n--- M ---");

    t_matrix M3 = matrixPower(&M, 3);
    printMatrix(&M3, "\n--- M^3 ---");

    t_matrix M7 = matrixPower(&M, 7);
    printMatrix(&M7, "\n--- M^7 ---");

    freeMatrix(&M);
    freeMatrix(&M3);
    freeMatrix(&M7);
}

int main() {
    char filename[100];

    printf("Nom du fichier de graphe : ");
    scanf("%s", filename);

    List_adj G = readGraph(filename);

    printf("\n--- Liste d'adjacence ---\n");
    DisplayListAdj(&G);

    printf("\n--- Verification Markov ---\n");
    Markov(G);

    printf("\n--- PARTIE 3 : tests matrices ---\n");
    run_part3_demo();

    return 0;
}
