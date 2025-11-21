#include <stdio.h>
#include "liste_adjacence.h"
#include "tarjan.h"

int main() {

    char filename[100];
    printf("Nom du fichier de graphe : ");
    scanf("%s", filename);

    List_adj G = readGraph(filename);

    printf("\n--- Liste d'adjacence ---\n");
    DisplayListAdj(&G);

    printf("\n--- Verification Markov ---\n");
    Markov(G);

    printf("\n--- Export Mermaid ---\n");
    ExportMermaid(G, "graphe_mermaid.txt");
    printf("Fichier genere : graphe_mermaid.txt\n");

    printf("\n--- Tarjan ---\n");
    t_partition P = tarjan(G);
    printPartition(P);
    return 0;
}
