#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "stack.h"

#define DEFAULT_SIZE 8

Stack stack_create(size_t element_size)
{
	Stack s;

	s.stack = (void*) malloc(DEFAULT_SIZE * element_size);
	s.logical_size = 0;
	s.physical_size = DEFAULT_SIZE;
	s.element_size = element_size;

	assert(s.stack);

	return s;
}

void stack_free(Stack* s)
{
	free(s->stack);

	s->logical_size = 0;
	s->physical_size = 0;
	s->element_size = 0;
}

void stack_push(Stack* s, void* elt_ptr)
{
	if(s->logical_size >= s->physical_size)
	{
		s->physical_size = 2*s->physical_size;
		s->stack = realloc(s->stack, s->physical_size * s->element_size);
	}

	memcpy(s->stack + (s->logical_size * s->element_size), elt_ptr, s->element_size);
	s->logical_size ++;
}

void stack_pull(Stack* s, void* elt_ptr)
{
	assert(!stack_isEmpty(s));

	if(elt_ptr)
		memcpy(elt_ptr, s->stack + (s->logical_size * s->element_size), s->element_size);
	s->logical_size --;
}

int stack_isEmpty(Stack* s)
{
	return s->logical_size <= 0;
}


size_t stack_size(Stack* s)
{
	return s->logical_size;
}