#ifndef LISTE_ADJACENCE_H
#define LISTE_ADJACENCE_H

#include <stdlib.h>
#include <stdio.h>

typedef struct t_cell {
    int sommet_d_arrive;
    float probabilite;
    struct t_cell *next;
} Cell;

typedef struct t_list {
    Cell *head;
} List;

typedef struct t_liste_d_adjacence {
    List *tab;
    int taille;
} List_adj;

// FONCTIONS
Cell* CreateCell(int sommet, float probabilite);
void AddCell(List* list, int sommet, float probabilite);
void DisplayList(List* list);
List_adj CreateEmptyListAdj(int taille);
void DisplayListAdj(List_adj* La);
void list_add_edge(List *L, int to, float p);
List_adj readGraph(const char *filename);
void Markov(List_adj adj);
char* getId(int num);
void ExportMermaid(List_adj G, const char *filename);

#endif