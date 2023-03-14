#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "bvhtree.h"
#include "vec3.h"
#include "bbox.h"
#include "material.h"
#include "light.h"
#include "typedef.h"

typedef struct Triangle Triangle;
struct Triangle
{
    long int vertices[3];
    long int uv[3];

    vec3 normal;
    vec3 dpdu;
    vec3 dpdv;

    long int material;
    long int group;
};

typedef struct SceneMetadata SceneMetadata;
struct SceneMetadata
{
    char **groups_name;
 
    size_t vertices_count;
    size_t vertices_tex_count;

    size_t triangles_count;
    size_t groups_count;
    
    size_t materials_count;
    size_t lights_count;

    BVHTree tree;
    
    bool use_uv;
};

typedef struct Scene Scene;
struct Scene
{
    vec3 *vertices;
    vec2 *vertices_tex;
    
    Triangle *triangles;

    Material *materials;
	Light* lights;
};

#endif