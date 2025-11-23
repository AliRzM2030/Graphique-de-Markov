#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tarjan.h"
#include "stack.h"
#include "hasse.h"

t_tarjan_vertex* initTarjanVertices(List_adj G) {
    t_tarjan_vertex* V = (t_tarjan_vertex*)malloc(G.taille * sizeof(t_tarjan_vertex));
    for (int i = 0; i < G.taille; i++) {
        V[i].id = i + 1;
        V[i].index = -1;
        V[i].lowlink = -1;
        V[i].onStack = 0;
    }
    return V;
}

void addVertexToClass(t_partition* part, int v) {
    t_classe* C = &part->classes[part->nb_classes - 1];
    C->vertices = (int*)realloc(C->vertices, sizeof(int) * (C->nb_vertices + 1));
    C->vertices[C->nb_vertices] = v;
    C->nb_vertices++;
}

void tarjanRecursive(int v, List_adj G, t_tarjan_vertex* V, Stack* P,
                     t_partition* part, int* index_tarjan) {

    V[v].index = *index_tarjan;
    V[v].lowlink = *index_tarjan;
    (*index_tarjan)++;

    push(P, v);
    V[v].onStack = 1;

    Cell* c = G.tab[v].head;
    while (c != NULL) {
        int w = c->sommet_d_arrive - 1;

        if (V[w].index == -1) {
            tarjanRecursive(w, G, V, P, part, index_tarjan);
            if (V[w].lowlink < V[v].lowlink)
                V[v].lowlink = V[w].lowlink;
        }
        else if (V[w].onStack && V[w].index < V[v].lowlink) {
            V[v].lowlink = V[w].index;
        }

        c = c->next;
    }

    if (V[v].lowlink == V[v].index) {
        part->classes = (t_classe*)realloc(part->classes,sizeof(t_classe) * (part->nb_classes + 1));
        t_classe* C = &part->classes[part->nb_classes];
        sprintf(C->name, "C%d", part->nb_classes + 1);
        C->vertices = NULL;
        C->nb_vertices = 0;
        part->nb_classes++;

        int w;
        do {
            w = pop(P);
            V[w].onStack = 0;
            addVertexToClass(part, w + 1);
        } while (w != v);
    }
}

t_partition tarjan(List_adj G) {
    t_partition part;
    part.classes = NULL;
    part.nb_classes = 0;

    t_tarjan_vertex* V = initTarjanVertices(G);
    Stack* P = createStack(G.taille);

    int index_tarjan = 0;

    for (int v = 0; v < G.taille; v++) {
        if (V[v].index == -1) {
            tarjanRecursive(v, G, V, P, &part, &index_tarjan);
        }
    }

    freeStack(P);
    free(V);
    return part;
}

void printPartition(t_partition part) {
    for (int i = 0; i < part.nb_classes; i++) {
        printf("%s : ", part.classes[i].name);
        for (int j = 0; j < part.classes[i].nb_vertices; j++) {
            printf("%d ", part.classes[i].vertices[j]);
        }
        printf("\n");
    }
}


//OPTIONNEL :
// ----------------------------------------------------------
// SUPPRESSION DES REDONDANCES DANS UNE CLASSE TARJAN
// ----------------------------------------------------------

#include <stdbool.h>

// Supprime les doublons dans une classe (C1, C2, ...)
// On garde l'ordre des sommets
void removeDuplicatesFromClass(t_classe* c) {
    if (c->nb_vertices <= 1) return;

    int* unique = malloc(c->nb_vertices * sizeof(int));
    int count = 0;

    for (int i = 0; i < c->nb_vertices; i++) {
        int v = c->vertices[i];
        bool exists = false;

        // Vérifier si déjà dans "unique"
        for (int j = 0; j < count; j++) {
            if (unique[j] == v) {
                exists = true;
                break;
            }
        }

        if (!exists) {
            unique[count++] = v;
        }
    }

    // Réécriture propre du tableau
    free(c->vertices);
    c->vertices = malloc(count * sizeof(int));
    for (int i = 0; i < count; i++) {
        c->vertices[i] = unique[i];
    }

    c->nb_vertices = count;
    free(unique);
}

// Supprime les redondances dans toutes les classes de la partition
void cleanPartition(t_partition* part) {
    for (int i = 0; i < part->nb_classes; i++) {
        removeDuplicatesFromClass(&(part->classes[i]));
    }
}
