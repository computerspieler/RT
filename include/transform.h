#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "typedef.h"

#include "matrix.h"

typedef struct Transform Transform;
struct Transform
{
    Matrix4x4 mat;
    Matrix4x4 matInv;
};

Transform transform_inverse(Transform t);
vec3 transform_apply_vector(vec3 p, Transform t);

#endif
