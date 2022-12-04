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
Transform transform_transpose(Transform t);
Transform transform_combine(Transform t1, Transform t2);

Transform transform_translate(vec3 delta);
Transform transform_scale(vec3 scale);

Transform transform_rotate_x(double angle);
Transform transform_rotate_y(double angle);
Transform transform_rotate_z(double angle);

vec3 transform_apply_vector(vec3 p, Transform t);
vec3 transform_apply_point(vec3 p, Transform t);

#endif
