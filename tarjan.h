//
// Created by melis on 22/11/2025.
//

#ifndef TARJAN_H
#define TARJAN_H
#include "liste_adjacence.h"
#include "stack.h"

typedef struct {
    int id;
    int index;
    int lowlink;
    int onStack;
} t_tarjan_vertex;

typedef struct {
    char name[10];
    int* vertices;
    int nb_vertices;
} t_classe;

typedef struct {
    t_classe* classes;
    int nb_classes;
} t_partition;



t_tarjan_vertex* initTarjanVertices(List_adj G);

void tarjanRecursive(int v, List_adj G, t_tarjan_vertex* V, Stack* P,
                     t_partition* part, int* index_tarjan);

t_partition tarjan(List_adj G);

void printPartition(t_partition part);

#endif