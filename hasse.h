//
// Created by melis on 23/11/2025.
//

#ifndef HASSE_H
#define HASSE_H

#include <stdlib.h>
#include <stdbool.h>
#include "tarjan.h"
#include "liste_adjacence.h"

typedef struct {
    int from;   // classe de départ
    int to;     // classe d'arrivée
} t_link;

typedef struct {
    t_link *data;
    int size;
    int capacity;
} t_link_array;

// Initialisation
void initLinkArray(t_link_array *arr, int capacity);

// Ajout de lien
void addLink(t_link_array *arr, int from, int to);

// Libération
void freeLinkArray(t_link_array *arr);

// Vérifie si un lien existe
bool linkExists(t_link_array *arr, int from, int to);

// Construit le tableau sommet → classe
int* buildVertexClassArray(t_partition part, int nb_vertices);

// Construit les liens entre classes
void buildClassLinks(List_adj G, t_partition part, t_link_array *links);

// Optionnel : suppression des liens transitifs
void removeTransitiveLinks(t_link_array *p_link_array);

// Export en format Mermaid
void exportHasseToMermaid(const char *filename, t_partition part, t_link_array *links);

#endif //HASSE_H