#ifndef _BBOX_H_
#define _BBOX_H_

#include "typedef.h"

typedef struct BBox3 BBox3;
struct BBox3
{
	vec3 min;
	vec3 max;
};

BBox3 bbox_b_union(BBox3 b1, BBox3 b2);
BBox3 bbox_p_union(BBox3 b, vec3 p);
int bbox_maximum_extent(BBox3 b);
vec3 bbox_p_offset(BBox3 b, vec3 p);
Float bbox_surface_area(BBox3 b);
vec3 bbox_center(BBox3 b);

BBox3 bbox_create(vec3 min, vec3 max);

#endif
