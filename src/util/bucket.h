#ifndef _BUFFER_H_
#define _BUFFER_H_

#include <sys/types.h>

typedef struct Bucket Bucket;
struct Bucket
{
	void *data;
	size_t length;
};

Bucket* bucket_create();
void* bucket_alloc(Bucket *b, size_t len);
void bucket_free(Bucket *b);

#endif