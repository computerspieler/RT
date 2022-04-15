#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "ray.h"

typedef struct Primitive Primitive;
struct Primitive
{
	void *metadata;
	double (*getCollidingTime) (Ray, void*);
};

#endif