#ifndef _RAY_H_
#define _RAY_H_

#define RAY_EPSILON 1e-3d

#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

struct Ray
{
	float3 o;				// Origine
	float3 d;				// Direction

	double t;				// Le temps
	double min_t, max_t;	// L'interval dans lequel le temps doit être compris
};

#endif