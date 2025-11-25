#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "liste_adjacence.h"
#include "tarjan.h"
#include "hasse.h"
#include "matrix.h"

void afficherMenu() {
    printf("\n===============================================================\n");
    printf("           MENU PRINCIPAL - GRAPHES DE MARKOV\n");
    printf("===============================================================\n");
    printf("1. Verifier si c'est un graphe de Markov\n");
    printf("2. Exporter le graphe au format Mermaid\n");
    printf("3. Algorithme de Tarjan (composantes fortement connexes)\n");
    printf("4. Diagramme de Hasse\n");
    printf("5. Analyser les caracteristiques du graphe\n");
    printf("6. Calculer la matrice et ses puissances\n");
    printf("7. Distribution stationnaire par classes\n");
    printf("8. Calculer la periode des classes\n");
    printf("9. Tout executer\n");
    printf("0. Quitter\n");
    printf("===============================================================\n");
    printf("Votre choix : ");
}

void analyserCaracteristiques(t_partition part, t_link_array links) {
    printf("\n===============================================================\n");
    printf("        ANALYSE DES CARACTERISTIQUES DU GRAPHE\n");
    printf("===============================================================\n");

    int *est_transitoire = calloc(part.nb_classes, sizeof(int));

    for (int i = 0; i < links.size; i++) {
        est_transitoire[links.data[i].from] = 1;
    }

    printf("\nCLASSES TRANSITOIRES :\n");
    int nb_transitoires = 0;
    for (int i = 0; i < part.nb_classes; i++) {
        if (est_transitoire[i]) {
            printf("  - %s : {", part.classes[i].name);
            for (int j = 0; j < part.classes[i].nb_vertices; j++) {
                printf("%d", part.classes[i].vertices[j]);
                if (j < part.classes[i].nb_vertices - 1) printf(", ");
            }
            printf("}\n");
            nb_transitoires++;
        }
    }
    if (nb_transitoires == 0) printf("  Aucune\n");

    printf("\nCLASSES PERSISTANTES :\n");
    for (int i = 0; i < part.nb_classes; i++) {
        if (!est_transitoire[i]) {
            printf("  - %s : {", part.classes[i].name);
            for (int j = 0; j < part.classes[i].nb_vertices; j++) {
                printf("%d", part.classes[i].vertices[j]);
                if (j < part.classes[i].nb_vertices - 1) printf(", ");
            }
            printf("}\n");
        }
    }

    printf("\nETATS ABSORBANTS :\n");
    int nb_absorbants = 0;
    for (int i = 0; i < part.nb_classes; i++) {
        if (!est_transitoire[i] && part.classes[i].nb_vertices == 1) {
            printf("  - Etat %d\n", part.classes[i].vertices[0]);
            nb_absorbants++;
        }
    }
    if (nb_absorbants == 0) printf("  Aucun\n");

    printf("\nIRREDUCTIBILITE :\n");
    if (part.nb_classes == 1) {
        printf("  [OK] Le graphe est IRREDUCTIBLE (une seule classe)\n");
    } else {
        printf("  [X] Le graphe n'est PAS irreductible (%d classes)\n", part.nb_classes);
    }

    printf("===============================================================\n");

    free(est_transitoire);
}

void calculerPuissancesMatrice(List_adj G) {
    printf("\n===============================================================\n");
    printf("          CALCUL DES PUISSANCES DE LA MATRICE\n");
    printf("===============================================================\n");

    t_matrix M = listAdjToMatrix(G);

    printf("\nMatrice de transition M :\n");
    printMatrix(&M, NULL);

    int choix;
    printf("\nQue souhaitez-vous calculer ?\n");
    printf("1. M^3\n");
    printf("2. M^7\n");
    printf("3. Convergence vers distribution stationnaire\n");
    printf("Votre choix : ");
    scanf("%d", &choix);

    if (choix == 1) {
        t_matrix M3 = matrixPower(&M, 3);
        printf("\nM^3 :\n");
        printMatrix(&M3, NULL);
        freeMatrix(&M3);
    } else if (choix == 2) {
        t_matrix M7 = matrixPower(&M, 7);
        printf("\nM^7 :\n");
        printMatrix(&M7, NULL);
        freeMatrix(&M7);
    } else if (choix == 3) {
        double epsilon = 0.01;
        int k = 1;
        t_matrix Mk = matrixPower(&M, k);
        t_matrix Mk_prev = createZeroMatrix(M.rows);

        printf("\nRecherche de convergence (epsilon = %.4f)...\n", epsilon);

        while (k < 1000) {
            copyMatrix(&Mk, &Mk_prev);
            freeMatrix(&Mk);
            k++;
            Mk = matrixPower(&M, k);

            double diff = diffMatrices(&Mk, &Mk_prev);

            if (diff < epsilon) {
                printf("\n[OK] Convergence atteinte a k = %d (diff = %.6f)\n", k, diff);
                printf("\nDistribution stationnaire M^%d :\n", k);
                printMatrix(&Mk, NULL);
                freeMatrix(&Mk);
                freeMatrix(&Mk_prev);
                freeMatrix(&M);
                return;
            }

            if (k % 100 == 0) {
                printf("  k = %d, diff = %.6f\n", k, diff);
            }
        }

        printf("\n[!] Pas de convergence detectee apres %d iterations\n", k);
        printf("Le graphe pourrait etre periodique ou non-convergent.\n");

        freeMatrix(&Mk);
        freeMatrix(&Mk_prev);
    }

    freeMatrix(&M);
}

void calculerDistributionsParClasses(List_adj G, t_partition part) {
    printf("\n===============================================================\n");
    printf("      DISTRIBUTIONS STATIONNAIRES PAR CLASSES\n");
    printf("===============================================================\n");

    t_matrix M = listAdjToMatrix(G);

    for (int c = 0; c < part.nb_classes; c++) {
        printf("\n--- Classe %s : {", part.classes[c].name);
        for (int i = 0; i < part.classes[c].nb_vertices; i++) {
            printf("%d", part.classes[c].vertices[i]);
            if (i < part.classes[c].nb_vertices - 1) printf(", ");
        }
        printf("} ---\n");

        t_matrix sub = subMatrix(M, part, c);

        if (sub.rows == 0) {
            printf("Erreur : sous-matrice invalide\n");
            continue;
        }

        printf("Sous-matrice :\n");
        printMatrix(&sub, NULL);

        double epsilon = 0.01;
        t_matrix Mk = matrixPower(&sub, 1);
        t_matrix Mk_prev = createZeroMatrix(sub.rows);

        int k;
        for (k = 2; k < 500; k++) {
            copyMatrix(&Mk, &Mk_prev);
            freeMatrix(&Mk);
            Mk = matrixPower(&sub, k);

            double diff = diffMatrices(&Mk, &Mk_prev);
            if (diff < epsilon) {
                printf("\n[OK] Distribution stationnaire (convergence a k=%d) :\n", k);
                printMatrix(&Mk, NULL);
                break;
            }
        }

        if (k >= 500) {
            printf("\n[!] Pas de convergence (classe possiblement periodique)\n");
        }

        freeMatrix(&Mk);
        freeMatrix(&Mk_prev);
        freeMatrix(&sub);
    }

    freeMatrix(&M);
}

void calculerPeriodes(List_adj G, t_partition part) {
    printf("\n===============================================================\n");
    printf("              CALCUL DES PERIODES\n");
    printf("===============================================================\n");

    t_matrix M = listAdjToMatrix(G);

    for (int c = 0; c < part.nb_classes; c++) {
        printf("\nClasse %s : {", part.classes[c].name);
        for (int i = 0; i < part.classes[c].nb_vertices; i++) {
            printf("%d", part.classes[c].vertices[i]);
            if (i < part.classes[c].nb_vertices - 1) printf(", ");
        }
        printf("}\n");

        t_matrix sub = subMatrix(M, part, c);

        if (sub.rows == 0) {
            printf("  Erreur : sous-matrice invalide\n");
            continue;
        }

        int periode = getPeriod(sub);
        printf("  Periode = %d\n", periode);

        if (periode == 1) {
            printf("  -> Classe aperiodique\n");
        } else {
            printf("  -> Classe periodique de periode %d\n", periode);
        }

        freeMatrix(&sub);
    }

    freeMatrix(&M);
}

void toutExecuter(List_adj G) {
    printf("\n===============================================================\n");
    printf("              EXECUTION COMPLETE\n");
    printf("===============================================================\n");

    printf("\n--- VERIFICATION MARKOV ---\n");
    Markov(G);

    ExportMermaid(G, "graphe_mermaid.mmd");
    printf("\n[OK] Graphe exporte vers 'graphe_mermaid.mmd'\n");

    printf("\n--- ALGORITHME DE TARJAN ---\n");
    t_partition part = tarjan(G);
    printPartition(part);

    printf("\n--- DIAGRAMME DE HASSE ---\n");
    t_link_array links;
    initLinkArray(&links, 10);
    buildClassLinks(G, part, &links);
    exportHasseToMermaid("hasse_mermaid.mmd", part, &links);
    printf("[OK] Diagramme de Hasse exporte vers 'hasse_mermaid.mmd'\n");

    analyserCaracteristiques(part, links);

    t_matrix M = listAdjToMatrix(G);
    printf("\n--- MATRICE DE TRANSITION ---\n");
    printMatrix(&M, NULL);

    calculerDistributionsParClasses(G, part);
    calculerPeriodes(G, part);

    freeMatrix(&M);
    freeLinkArray(&links);
    for (int i = 0; i < part.nb_classes; i++) {
        free(part.classes[i].vertices);
    }
    free(part.classes);
}

int main() {
    char filename[256];

    printf("===============================================================\n");
    printf("        PROJET GRAPHES DE MARKOV - TI301\n");
    printf("===============================================================\n\n");

    printf("Entrez le nom du fichier graphe : ");
    scanf("%s", filename);

    List_adj G = readGraph(filename);

    printf("\n[OK] Graphe charge : %d sommets\n", G.taille);

    printf("\n=== DEBUG : Contenu du graphe ===\n");
    for (int i = 0; i < G.taille; i++) {
        printf("Sommet %d : ", i+1);
        Cell *c = G.tab[i].head;
        while (c != NULL) {
            printf("-> %d (%.2f) ", c->sommet_d_arrive, c->probabilite);
            c = c->next;
        }
        printf("\n");
    }
    printf("=================================\n");

    int choix;
    t_partition part = {NULL, 0};
    t_link_array links = {NULL, 0, 0};
    int tarjan_execute = 0;
    do {
    afficherMenu();
    scanf("%d", &choix);

    switch(choix) {
        case 1:
            printf("\n--- VERIFICATION MARKOV ---\n");
            Markov(G);
            break;

        case 2:
            ExportMermaid(G, "graphe_mermaid.mmd");
            printf("\n[OK] Graphe exporte vers 'graphe_mermaid.mmd'\n");
            break;

        case 3:
            if (tarjan_execute) {
                for (int i = 0; i < part.nb_classes; i++) {
                    free(part.classes[i].vertices);
                }
                free(part.classes);
            }
            printf("\n--- ALGORITHME DE TARJAN ---\n");
            part = tarjan(G);
            tarjan_execute = 1;
            printPartition(part);
            break;

        case 4:
            if (!tarjan_execute) {
                printf("\n[!] Veuillez d'abord executer Tarjan (option 3)\n");
                break;
            }
            if (links.data) freeLinkArray(&links);
            initLinkArray(&links, 10);
            buildClassLinks(G, part, &links);
            exportHasseToMermaid("hasse_mermaid.mmd", part, &links);
            printf("\n[OK] Diagramme de Hasse exporte vers 'hasse_mermaid.mmd'\n");
            break;

        case 5:
            if (!tarjan_execute) {
                printf("\n[!] Veuillez d'abord executer Tarjan (option 3)\n");
                break;
            }
            if (!links.data || links.size == 0) {
                initLinkArray(&links, 10);
                buildClassLinks(G, part, &links);
            }
            analyserCaracteristiques(part, links);
            break;

        case 6:
            calculerPuissancesMatrice(G);
            break;

        case 7:
            if (!tarjan_execute) {
                printf("\n[!] Veuillez d'abord executer Tarjan (option 3)\n");
                break;
            }
            calculerDistributionsParClasses(G, part);
            break;

        case 8:
            if (!tarjan_execute) {
                printf("\n[!] Veuillez d'abord executer Tarjan (option 3)\n");
                break;
            }
            calculerPeriodes(G, part);
            break;

        case 9:
            if (tarjan_execute) {
                for (int i = 0; i < part.nb_classes; i++) {
                    free(part.classes[i].vertices);
                }
                free(part.classes);
                if (links.data) freeLinkArray(&links);
            }
            toutExecuter(G);
            tarjan_execute = 0;
            part.classes = NULL;
            part.nb_classes = 0;
            links.data = NULL;
            break;

        case 0:
            printf("\nAu revoir !\n");
            break;

        default:
            printf("\n[!] Choix invalide\n");
    }
} while (choix != 0);

if (tarjan_execute) {
    for (int i = 0; i < part.nb_classes; i++) {
        free(part.classes[i].vertices);
    }
    free(part.classes);
}
if (links.data) freeLinkArray(&links);

for (int i = 0; i < G.taille; i++) {
    Cell *c = G.tab[i].head;
    while (c) {
        Cell *tmp = c;
        c = c->next;
        free(tmp);
    }
}
free(G.tab);

return 0;
}