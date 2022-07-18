#ifndef _TREE_H_
#define _TREE_H_

#include <stdbool.h>

typedef struct Node Node;
typedef struct Edge Edge;

struct Edge
{
    char c;
    Node *son;
};

struct Node
{
    bool end;
    int value;
    Edge* edges;
    int edges_count;
};

void add_word(Node* n, char *text, int value);
Node* get_node(Node* n, char *text);
void print_tree(Node *n, int indent);
void free_tree(Node *n);

#endif