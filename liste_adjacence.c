#include "liste_adjacence.h"

// -------------------------------
// Création d’une cellule
// -------------------------------
Cell* CreateCell(int sommet , float probabilite){
    Cell *newcell = (Cell*)malloc(sizeof(Cell));
    if (newcell == NULL) {
        perror("Erreur d'allocation mémoire");
        exit(EXIT_FAILURE);
    }
    newcell->sommet_d_arrive = sommet;
    newcell->probabilite = probabilite;
    newcell->next = NULL;
    return newcell;
}

// -------------------------------
// Ajout d'une cellule dans une liste
// -------------------------------
void AddCell(List* list , int sommet , float probabilite){
    Cell* newCell = CreateCell(sommet, probabilite);
    newCell->next = list->head;
    list->head = newCell;
}

// -------------------------------
// Affichage d'une liste
// -------------------------------
void DisplayList(List* list){
    Cell* curr = list->head;
    while(curr != NULL){
        printf("(%d, %.2f) ", curr->sommet_d_arrive, curr->probabilite);
        curr = curr->next;
    }
}

// -------------------------------
// Liste d'adjacence vide
// -------------------------------
List_adj CreateEmptyListAdj(int taille) {
    List_adj la;
    la.taille = taille;
    la.tab = (List*)malloc(taille * sizeof(List));

    if (la.tab == NULL) {
        perror("Erreur allocation du tableau");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < taille; ++i) {
        la.tab[i].head = NULL;
    }

    return la;
}

// -------------------------------
// Affichage liste d'adjacence
// -------------------------------
void DisplayListAdj(List_adj* La){
    for(int i = 0 ; i < La->taille ; i++){
        printf("Sommet %d : ", i+1);
        DisplayList(&La->tab[i]);  // IMPORTANT : passer l’adresse
        printf("\n");
    }
}

// -------------------------------
// Ajout d'une arête
// -------------------------------
void list_add_edge(List *L, int to, float p) {

    if (p <= 0.0f || p > 1.0f) {
        fprintf(stderr, "Probabilité invalide : %.2f\n", p);
        exit(EXIT_FAILURE);
    }

    Cell *newCell = malloc(sizeof(Cell));
    if (!newCell){
        perror("Erreur d'allocation");
        exit(EXIT_FAILURE);
    }

    newCell->sommet_d_arrive = to;
    newCell->probabilite = p;

    newCell->next = L->head;
    L->head = newCell;
}

// -------------------------------
// Lecture d'un fichier
// -------------------------------
List_adj readGraph(const char *filename) {
    FILE *file = fopen(filename, "rt");
    if (!file) {
        perror("Could not open file for reading");
        exit(EXIT_FAILURE);
    }

    int nbvert, depart, arrivee;
    float proba;

    fscanf(file, "%d", &nbvert);
    List_adj G = CreateEmptyListAdj(nbvert);

    while (fscanf(file, "%d %d %f", &depart, &arrivee, &proba) == 3) {

        if (depart < 1 || depart > nbvert || arrivee < 1 || arrivee > nbvert) {
            fprintf(stderr, "Sommet hors bornes: %d -> %d (n=%d)\n",
                    depart, arrivee, nbvert);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        list_add_edge(&G.tab[depart - 1], arrivee, proba);
    }

    fclose(file);
    return G;
}

// -------------------------------
// Vérification Markov
// -------------------------------
void Markov(List_adj adj) {
    int ok = 1;

    for (int i = 0; i < adj.taille; i++) {
        float sum = 0;
        Cell *m = adj.tab[i].head;

        while (m != NULL) {
            sum += m->probabilite;
            m = m->next;
        }

        if (sum < 0.99f || sum > 1.001f) {
            printf("Sommet %d NON valide : somme = %.2f\n", i+1, sum);
            ok = 0;
        }
    }

    if (ok) printf("Le graphe est un graphe de Markov\n");
}

// -------------------------------
// getId : A, B, C, ... AA, AB...
// -------------------------------
char* getId(int num) {
    static char buf[10];
    int index = 0;
    int n = num;

    while (n > 0) {
        n--;
        buf[index++] = 'A' + (n % 26);
        n /= 26;
    }

    buf[index] = '\0';

    for (int i = 0; i < index/2; i++) {
        char tmp = buf[i];
        buf[i] = buf[index - i - 1];
        buf[index - i - 1] = tmp;
    }

    return buf;
}

// -------------------------------
// Exportation version Mermaid
// -------------------------------
void ExportMermaid(List_adj G, const char *filename) {

    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Erreur création du fichier Mermaid");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "---\nconfig:\n  layout: elk\n  theme: neo\n  look: neo\n---\n");
    fprintf(f, "flowchart LR\n");

    for (int i = 0; i < G.taille; i++){
        fprintf(f, "%s((%d))\n", getId(i+1), i+1);
    }

    for (int i = 0; i < G.taille; i++){
        Cell *c = G.tab[i].head;
        while(c){
            fprintf(f, "%s -->|%.2f| %s\n",
                getId(i+1),
                c->probabilite,
                getId(c->sommet_d_arrive));
            c = c->next;
        }
    }

    fclose(f);
}