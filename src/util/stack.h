#ifndef _STACK_H_
#define _STACK_H_

typedef struct {
	void *stack;
	
	size_t physical_size;
	size_t logical_size;
	size_t element_size;
} Stack;

Stack stack_create(size_t element_size);
void stack_free(Stack*);
void stack_push(Stack*, void*);
void stack_pull(Stack*, void*);
size_t stack_size(Stack*);
int stack_isEmpty(Stack*);

#endif
