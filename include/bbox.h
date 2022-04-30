#ifndef _BBOX_H_
#define _BBOX_H_

#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

typedef struct BBox3 BBox3;
struct BBox3
{
	float3 min;
	float3 max;
};

#endif