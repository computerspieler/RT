#ifndef _FORM_H_
#define _FORM_H_

#include "vec3.h"
#include "primitive.h"

typedef struct SphereMetadata SphereMetadata;
struct SphereMetadata
{
	Vec3 pos;
	double radius;
};

Primitive *createSphere(SphereMetadata);

#endif