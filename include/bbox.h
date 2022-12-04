#ifndef _BBOX_H_
#define _BBOX_H_

#include "typedef.h"

typedef struct BBox3 BBox3;
struct BBox3
{
	vec3 min;
	vec3 max;
};

BBox3 bbox_create(vec3 min, vec3 max);
vec3 bbox_corner(BBox3 b, int corner);
vec3 bbox_center(BBox3 b);
BBox3 bbox_p_union(BBox3 b, vec3 p);
BBox3 bbox_b_union(BBox3 b1, BBox3 b2);
BBox3 bbox_b_intersection(BBox3 b1, BBox3 b2);
bool bbox_b_overlap(BBox3 b1, BBox3 b2);
bool bbox_p_inside(vec3 p, BBox3 b);
bool bbox_p_insideExclusive(vec3 p, BBox3 b);
BBox3 bbox_p_expand(BBox3 b, double delta);
vec3 bbox_diagonal(BBox3 b);
double bbox_surface_area(BBox3 b);
double bbox_surface_volume(BBox3 b);
int bbox_maximum_extent(BBox3 b);
vec3 bbox_lerp(BBox3 b, vec3 p);
vec3 bbox_p_offset(BBox3 b, vec3 p);

#endif
