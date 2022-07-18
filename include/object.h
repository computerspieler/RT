#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "double3.h"
#include "bbox.h"
#include "material.h"


#include "typedef.h"

typedef struct Triangle Triangle;
struct Triangle
{
    uint vertices[3];
    uint uv[3];
    uint normals[3];

    uint material;
    uint group;
};

typedef struct ObjectMetadata ObjectMetadata;
struct ObjectMetadata
{
    char **groups_name;
 
    size_t vertices_count;
    size_t vertices_normal_count;
    size_t vertices_tex_count;

    size_t triangles_count;
    size_t groups_count;
    
    size_t materials_count;

    BBox3 bounds;
};

typedef struct Object Object;
struct Object
{
    double3 *vertices;
    double3 *vertices_normal;
    double2 *vertices_tex;
    
    Triangle *triangles;

    Material *materials;
};

#endif