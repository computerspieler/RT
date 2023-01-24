#ifndef _RAY_H_
#define _RAY_H_

#include "material.h"
#include "scene.h"
#include "typedef.h"

typedef struct Ray Ray;
struct Ray
{
	vec3 origin;				// Origine
	vec3 direction;			// Direction

	double min_t, max_t;	// L'interval dans lequel le temps doit être compris
};

typedef struct SurfaceInteraction SurfaceInteraction;
struct SurfaceInteraction
{
	vec3 p;		// Point
	vec3 n;		// Normal
	double2 uv;

	vec3 dpdu;
	vec3 dpdv;

	double time;
	int triangle_id;
	Material m;
};


#endif
