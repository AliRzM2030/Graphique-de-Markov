#include <stdio.h>
#include <string.h>
#include "hasse.h"

void initLinkArray(t_link_array *arr, int capacity) {
    arr->data = (t_link*)malloc(capacity * sizeof(t_link));
    arr->size = 0;
    arr->capacity = capacity;
}

void addLink(t_link_array *arr, int from, int to) {
    if (arr->size >= arr->capacity) {
        int newCap = (arr->capacity == 0) ? 4 : arr->capacity * 2;
        arr->data = (t_link*)realloc(arr->data, newCap * sizeof(t_link));
        arr->capacity = newCap;
    }
    arr->data[arr->size].from = from;
    arr->data[arr->size].to = to;
    arr->size++;
}

void freeLinkArray(t_link_array *arr) {
    free(arr->data);
    arr->data = NULL;
    arr->size = 0;
    arr->capacity = 0;
}

bool linkExists(t_link_array *arr, int from, int to) {
    for (int i = 0; i < arr->size; i++) {
        if (arr->data[i].from == from && arr->data[i].to == to)
            return true;
    }
    return false;
}

int* buildVertexClassArray(t_partition part, int nb_vertices) {
    int *class_of = (int*)malloc(nb_vertices * sizeof(int));
    
    for (int c = 0; c < part.nb_classes; c++) {
        for (int j = 0; j < part.classes[c].nb_vertices; j++) {
            int v = part.classes[c].vertices[j] - 1; // CORRECTION: sommets base 1
            if(v >= 0 && v < nb_vertices) {
                class_of[v] = c;
            }
        }
    }
    return class_of;
}
void buildClassLinks(List_adj G, t_partition part, t_link_array *links) {
    int *class_of = buildVertexClassArray(part, G.taille);
    for (int i = 0; i < G.taille; i++) {
        int Ci = class_of[i];
        Cell *c = G.tab[i].head;

        while (c != NULL) {
            int j = c->sommet_d_arrive - 1;  // CORRECTION: base 1 → base 0

            if (j >= 0 && j < G.taille) {
                int Cj = class_of[j];

                if (Ci != Cj && !linkExists(links, Ci, Cj)) {
                    addLink(links, Ci, Cj);
                }
            }
            c = c->next;
        }
    }

    free(class_of);
}
/* CORRECTION: Export Mermaid pour Hasse */
void exportHasseToMermaid(const char *filename, t_partition part, t_link_array *links) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Impossible de créer le fichier Mermaid");
        return;
    }
    // En-tête Mermaid
    fprintf(f, "---\n");
    fprintf(f, "config:\n");
    fprintf(f, "  layout: elk\n");
    fprintf(f, "  theme: neo\n");
    fprintf(f, "  look: neo\n");
    fprintf(f, "---\n");
    fprintf(f, "graph TD\n");

    // Définition des nœuds (classes)
    for (int i = 0; i < part.nb_classes; i++) {
        fprintf(f, "    %s[\"{", part.classes[i].name);
        for (int j = 0; j < part.classes[i].nb_vertices; j++) {
            fprintf(f, "%d", part.classes[i].vertices[j]);
            if (j < part.classes[i].nb_vertices - 1)
                fprintf(f, ",");
        }
        fprintf(f, "}\"]\n");
    }

    // Liens entre classes
    for (int i = 0; i < links->size; i++) {
        fprintf(f, "    C%d --> C%d\n",
                links->data[i].from + 1,
                links->data[i].to + 1);
    }

    fclose(f);
}