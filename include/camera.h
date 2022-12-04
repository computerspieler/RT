#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "typedef.h"

#include "transform.h"

typedef struct Camera Camera;
struct Camera
{
    vec3 pos;
    vec3 rot;

    double near;
    double max_t;
    double fov;
    int2 viewport;

    Transform rotation_transform;
};

#endif
