#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "typedef.h"

#include "matrix.h"

typedef struct Camera Camera;
struct Camera
{
    double3 pos;
    double3 rot;

    double near;
    double max_t;
    double fov;
    int2 viewport;

    Matrix4x4 rotation_matrix;
};

#endif
