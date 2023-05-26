#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "typedef.h"

typedef enum MaterialType MaterialType;
enum MaterialType
{
    MATERIAL_NONE = 0,
    MATERIAL_LAMBERTIAN = 1,
    MATERIAL_MIRROR = 3,
    MATERIAL_GLASS = 4,
};

typedef struct MaterialLambertian MaterialLambertian;
struct MaterialLambertian
{
    Float R;
};

typedef struct MaterialOren_Nayar MaterialOren_Nayar;
struct MaterialOren_Nayar
{
    Float sigma;
    Float R;
};

typedef struct MaterialGlass MaterialGlass;
struct MaterialGlass
{
    Float IOR;
};

typedef struct MaterialComplex MaterialComplex;
struct MaterialComplex
{
    Float specular_roughness;
    Float specular_IOR;
    Float metalness;
    Float transmission;
    Float transmission_dispersion;
};

typedef struct Material Material;
struct Material
{
    MaterialType type;

    union {
        MaterialLambertian l;
        MaterialOren_Nayar on;
        MaterialGlass g;
        MaterialComplex c;
    };
};

#endif
