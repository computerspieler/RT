#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "typedef.h"

typedef struct Material Material;
struct Material
{
    double specular;
    double density;
    double transparency;
    
    double3 ambient_color;
    double3 specular_color;
    double3 emissive_color;
    double3 diffuse_color;
};

#endif
