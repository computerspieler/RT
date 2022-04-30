#ifndef _SCENE_H_
#define _SCENE_H_

#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

typedef struct Material Material;
struct Material
{
	float3 color;
};

typedef struct Sphere Sphere;
struct Sphere
{
	float3 center;
	float radius;
	int material;
};

typedef struct Scene Scene;
struct Scene {
	int spheres_count;
	int materials_count;

	Sphere *spheres;
	Material *materials;
};

#endif