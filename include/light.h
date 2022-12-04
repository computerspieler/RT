#ifndef _LIGHT_H_
#define _LIGHT_H_

#include "typedef.h"
#include "spectrum.h"

typedef enum LightType LightType;
enum LightType
{
    LIGHT_POINT
};

typedef struct Light Light;
struct Light
{
    vec3 pos;
    LightType type;
    double I;
    Spectrum light;
};


#endif
