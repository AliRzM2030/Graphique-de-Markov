//
// Created by melis on 22/11/2025.
//


#ifndef STACK_H
#define STACK_H

typedef struct {
    int* data;
    int top;
    int capacity;
} Stack;

Stack* createStack(int capacity);
void push(Stack* s, int value);
int pop(Stack* s);
int peek(Stack* s);
int isEmpty(Stack* s);
void freeStack(Stack* s);

#endif