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

Transform transform_translate(double3 delta);
Transform transform_scale(double3 scale);

Transform transform_rotate_x(double angle);
Transform transform_rotate_y(double angle);
Transform transform_rotate_z(double angle);

double3 transform_apply_vector(double3 p, Transform t);
double3 transform_apply_point(double3 p, Transform t);

#endif
