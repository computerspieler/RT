#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "typedef.h"
#include "spectrum.h"
#include "map.h"

typedef enum MaterialType MaterialType;
enum MaterialType
{
    MATERIAL_NONE = 0,
    MATERIAL_LAMBERTIAN = 1,
    MATERIAL_OREN_NAYAR = 2,
    MATERIAL_MIRROR = 3,
    MATERIAL_GLASS = 4,
};

typedef struct MaterialLambertian MaterialLambertian;
struct MaterialLambertian
{
    double rho;
};

typedef struct MaterialOren_Nayar MaterialOren_Nayar;
struct MaterialOren_Nayar
{
    double sigma;
    double R;
};

typedef struct MaterialGlass MaterialGlass;
struct MaterialGlass
{
    double IOR;
};

typedef struct MaterialComplex MaterialComplex;
struct MaterialComplex
{
    vec3 specular_color;
    double specular_roughness;
    double specular_IOR;
    double metalness;
    double transmission;
    vec3 transmission_color;
    double transmission_dispersion;
};

typedef struct Material Material;
struct Material
{
    MaterialType type;
	double Le;
    
	Map normal_map;

    union {
        MaterialLambertian l;
        MaterialOren_Nayar on;
        MaterialGlass g;
        MaterialComplex c;
    };
};

#endif
