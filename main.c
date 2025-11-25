#include <stdio.h>
#include <stdlib.h>
#include "liste_adjacence.h"
#include "tarjan.h"
#include "hasse.h"
#include "matrix.h"

void freeAdjList(t_adj_list *A, int n) {
    for (int i = 0; i < n; i++)
        free(A[i].edges);
    free(A);
}

int main() {
    List_adj G;
    t_partition part;
    t_link_array links;
    char filename[256];

    printf("Entrez le nom du fichier du graphe : ");
    scanf("%255s", filename);

    G = readGraph(filename);

    int choix;

    do {
        printf("\n===========================================\n");
        printf("              MENU PRINCIPAL\n");
        printf("===========================================\n");
        printf("1. Verifier si c'est un graphe de Markov\n");
        printf("2. Exporter le graphe au format Mermaid\n");
        printf("3. Algorithme de Tarjan\n");
        printf("4. Diagramme de Hasse\n");
        printf("5. Caracteristiques du graphe\n");
        printf("6. Matrice stationnaire\n");
        printf("7. Distribution stationnaire par classes\n");
        printf("8. Tout executer\n");
        printf("0. Quitter\n");
        printf("===========================================\n");
        printf("Votre choix : ");
        scanf("%d", &choix);

        switch (choix) {

        case 1:
            Markov(G);
            break;

        case 2:
            ExportMermaid(G, "graph_mermaid.md");
            printf("Fichier graph_mermaid.md généré.\n");
            break;

        case 3:
            part = tarjan(G);
            printPartition(part);
            break;

        case 4:
            part = tarjan(G);
            initLinkArray(&links, 10);
            buildClassLinks(G, part, &links);
            removeTransitiveLinks(&links);
            exportHasseToMermaid("hasse.md", part, &links);
            freeLinkArray(&links);
            printf("Fichier hasse.md généré.\n");
            break;

        case 5: {
            part = tarjan(G);
            initLinkArray(&links, 10);
            buildClassLinks(G, part, &links);

            printf("\nCaractéristiques :\n");

            for (int c = 0; c < part.nb_classes; c++) {
                int outgoing = 0;

                for (int i = 0; i < links.size; i++)
                    if (links.data[i].from == c)
                        outgoing = 1;

                printf("Classe C%d : %s\n", c+1,
                        outgoing ? "transitoire" : "persistante");

                if (!outgoing && part.classes[c].nb_vertices == 1)
                    printf("   -> Etat absorbant : %d\n", part.classes[c].vertices[0]);
            }

            printf("Graphe %s irreductible\n",
                (part.nb_classes == 1 ? "" : "non"));

            freeLinkArray(&links);
        } break;

        case 6: {
            t_adj_list *A = convertListAdj_to_adjList(G);
            t_matrix M = adjacencyListToMatrix(A, G.taille);

            t_matrix Mk = matrixPower(&M, 1);
            double diff;
            int k = 1;

            do {
                t_matrix next = matrixPower(&M, k+1);
                diff = diffMatrices(&Mk, &next);
                freeMatrix(&Mk);
                Mk = next;
                k++;
            } while (diff > 0.01);

            printMatrix(&Mk, "Matrice stationnaire :");

            freeMatrix(&Mk);
            freeMatrix(&M);
            freeAdjList(A, G.taille);
        } break;

        case 7: {
            part = tarjan(G);

            t_adj_list *A = convertListAdj_to_adjList(G);
            t_matrix M = adjacencyListToMatrix(A, G.taille);

            for (int c = 0; c < part.nb_classes; c++) {
                printf("\n--- Classe C%d ---\n", c+1);

                t_matrix sub = subMatrix(M, part, c);
                int per = getPeriod(sub);
                printf("Période : %d\n", per);

                t_matrix Mk = matrixPower(&sub, 20);
                printMatrix(&Mk, "Stationnaire approximative :");

                freeMatrix(&sub);
                freeMatrix(&Mk);
            }

            freeMatrix(&M);
            freeAdjList(A, G.taille);
        } break;

        case 8: {
    printf("\n===== TOUT EXECUTER =====\n");

    /* [1] Vérifier Markov */
    printf("\n[1] Vérification Markov\n");
    Markov(G);

    /* [2] Export Mermaid du graphe */
    printf("\n[2] Export Mermaid du graphe\n");
    ExportMermaid(G, "graph_mermaid.md");
    printf("Fichier graph_mermaid.md généré.\n");

    /* [3] Tarjan : classes fortement connexes */
    printf("\n[3] Algorithme de Tarjan\n");
    part = tarjan(G);
    printPartition(part);

    /* [4] Diagramme de Hasse */
    printf("\n[4] Diagramme de Hasse\n");
    initLinkArray(&links, 10);
    buildClassLinks(G, part, &links);
    removeTransitiveLinks(&links);
    exportHasseToMermaid("hasse.md", part, &links);
    printf("Fichier hasse.md généré.\n");
    freeLinkArray(&links);

    /* [5] Caractéristiques du graphe */
    printf("\n[5] Caractéristiques du graphe\n");

    initLinkArray(&links, 10);
    buildClassLinks(G, part, &links);

    int is_irreducible = (part.nb_classes == 1);

    for (int c = 0; c < part.nb_classes; c++) {
        int outgoing = 0;

        for (int i = 0; i < links.size; i++) {
            if (links.data[i].from == c)
                outgoing = 1;
        }

        printf("Classe C%d : %s\n",
               c + 1,
               outgoing ? "transitoire" : "persistante");

        if (!outgoing && part.classes[c].nb_vertices == 1) {
            printf("   -> Etat absorbant : %d\n",
                   part.classes[c].vertices[0]);
        }
    }

    printf("Le graphe est %sirréductible\n",
           is_irreducible ? "" : "non ");

    freeLinkArray(&links);

    /* [6] Matrice stationnaire globale */
    printf("\n[6] Matrice stationnaire (M^k)\n");

    t_adj_list *A = convertListAdj_to_adjList(G);
    t_matrix M = adjacencyListToMatrix(A, G.taille);

    t_matrix Mk = matrixPower(&M, 1);
    double diff;
    int k = 1;

    do {
        t_matrix next = matrixPower(&M, k + 1);
        diff = diffMatrices(&Mk, &next);
        freeMatrix(&Mk);
        Mk = next;
        k++;
    } while (diff > 0.01);

    printMatrix(&Mk, "Matrice stationnaire approximative :");
    freeMatrix(&Mk);

    /* [7] Distributions stationnaires par classes */
    printf("\n[7] Distributions stationnaires par classes\n");

    for (int c = 0; c < part.nb_classes; c++) {
        printf("\n--- Classe C%d ---\n", c + 1);

        t_matrix sub = subMatrix(M, part, c);
        int per = getPeriod(sub);
        printf("Période : %d\n", per);

        t_matrix Mk_sub = matrixPower(&sub, 20);
        printMatrix(&Mk_sub, "Distribution stationnaire approx :");

        freeMatrix(&sub);
        freeMatrix(&Mk_sub);
    }

    freeMatrix(&M);
    freeAdjList(A, G.taille);

    printf("\n===== FIN DE TOUT EXECUTER =====\n");
} break;

        case 0:
            printf("Au revoir !\n");
            break;

        default:
            printf("Choix invalide !\n");
        }

    } while (choix != 0);

    return 0;
}