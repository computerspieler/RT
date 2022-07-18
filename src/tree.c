#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tree.h"

void add_word(Node* n, char *text, int value)
{
    int i;

    if(text[0] == 0) {
        n->value = value;
        n->end = true;
        return;
    }

    for(i = 0; i < n->edges_count; i++)
        if(text[0] == n->edges[i].c) {
            add_word(n->edges[i].son, text + 1, value);
            return;
        }
    
    n->edges_count ++;
    n->edges = realloc(n->edges, sizeof(Edge) * n->edges_count);
    n->edges[n->edges_count - 1].c = text[0];
    n->edges[n->edges_count - 1].son = malloc(sizeof(Node));
    bzero(n->edges[n->edges_count - 1].son, sizeof(Node));

    add_word(n->edges[n->edges_count - 1].son, text + 1, value);
}

Node* get_node(Node* n, char *text)
{
    int i;

    if(text[0] == 0)
        return n;

    for(i = 0; i < n->edges_count; i++)
        if(text[0] == n->edges[i].c)
            return get_node(n->edges[i].son, text + 1);

    return NULL;
}
 
void free_tree(Node *n)
{
    int i;

    for(i = 0; i < n->edges_count; i++)
        free_tree(n->edges[i].son);
    
    free(n);
}

void print_tree(Node *n, int indent)
{
    int i;

    if(n->edges_count != 1) {
        putchar(n->end ? '*' : '-');
        putchar('\n');
        for(i = 0; i < n->edges_count; i++) {
            for(int j = 0; j < indent; j ++)
                putchar(' ');
            putchar(n->edges[i].c);
            print_tree(n->edges[i].son, indent + 1);
        }
    } else {
        putchar(n->end ? '*' : '\0');
        putchar(n->edges[0].c);
        print_tree(n->edges[0].son, indent + 1);
    }

}