//
// Created by melis on 22/11/2025.
//


#include <stdlib.h>
#include <stdio.h>
#include "stack.h"

Stack* createStack(int capacity) {
    Stack* s = (Stack*)malloc(sizeof(Stack));
    if (!s) {
        perror("stack allocation error");
        exit(EXIT_FAILURE);
    }

    s->data = (int*)malloc(sizeof(int) * capacity);
    if (!s->data) {
        perror("stack data allocation error");
        exit(EXIT_FAILURE);
    }

    s->top = -1;
    s->capacity = capacity;

    return s;
}

void push(Stack* s, int value) {
    if (s->top + 1 >= s->capacity) {
        s->capacity *= 2;
        s->data = (int*)realloc(s->data, sizeof(int) * s->capacity);
        if (!s->data) {
            perror("stack realloc error");
            exit(EXIT_FAILURE);
        }
    }

    s->data[++s->top] = value;
}

int pop(Stack* s) {
    if (s->top < 0) {
        printf("stack underflow\n");
        return -1;
    }
    return s->data[s->top--];
}

int peek(Stack* s) {
    if (s->top < 0) {
        return -1;
    }
    return s->data[s->top];
}

int isEmpty(Stack* s) {
    return (s->top < 0);
}

void freeStack(Stack* s) {
    free(s->data);
    free(s);
}