#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "bucket.h"

Bucket* bucket_create()
{
	Bucket* b = malloc(sizeof(Bucket));
	if(!b) {
		perror("malloc");
		return NULL;
	}

	memset(b, 0, sizeof(Bucket));

	return b;
}

void* bucket_alloc(Bucket *b, size_t len)
{
	void *output;

	if(!b)
		return NULL;

	if(b->length) {
		output = (b->data + b->length);
		b->length += len;
		b->data = realloc(b->data, b->length);
	} else {
		b->length = len;
		output = b->data = malloc(len);
	}

	return output;
}

void bucket_free(Bucket *b)
{
	if(!b) return;

	free(b->data);
	free(b);
}