#ifndef _LIGHT_H_
#define _LIGHT_H_

#include "typedef.h"

typedef enum LightType LightType;
enum LightType
{
    LIGHT_POINT,
	LIGHT_AREA
};

typedef struct Light Light;
struct Light
{
    vec3 pos;
	vec3 area_n;
	vec3 area_tan;
	vec3 area_bitan;
	Float area_width;
	Float area_height;
    LightType type;
    Float I;
};

#endif
