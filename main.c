#include <stdio.h>
#include "liste_adjacence.h"
#include "tarjan.h"
#include "hasse.h"
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

    t_partition Pa = tarjan(G);
    printPartition(Pa);

    // Construire les liens entre classes
    t_link_array links;
    initLinkArray(&links, 8);

    buildClassLinks(G, Pa, &links);

    // Optionnel :
    removeTransitiveLinks(&links);

    // Export Hasse
    exportHasseToMermaid("hasse_mermaid.txt", Pa, &links);

    freeLinkArray(&links);
}
