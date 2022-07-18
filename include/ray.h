#ifndef _RAY_H_
#define _RAY_H_

#define RAY_EPSILON 1e-4

#include "typedef.h"

typedef struct Ray Ray;
struct Ray
{
	double3 origin;				// Origine
	double3 direction;			// Direction

	double t;				// Le temps
	double min_t, max_t;	// L'interval dans lequel le temps doit être compris
};

typedef struct SurfaceInteraction SurfaceInteraction;
struct SurfaceInteraction
{
	double3 p;		// Point
	double3 n;		// Normal
	double3 dpdu;	// Derivative of p on u
	double3 dpdv;	// Derivative of p on v
	double3 dndu;
	double3 dndv;
	double time;
};

typedef struct RayDifferential RayDifferential;
struct RayDifferential
{
	Ray ray;

	bool has_differential;
	double3 rx_origin, ry_origin;
	double3 rx_direction, ry_direction;
};

RayDifferential ray_scale_differentials(RayDifferential r, double s);

#endif
