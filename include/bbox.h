#ifndef _BBOX_H_
#define _BBOX_H_

#include "typedef.h"

typedef struct BBox3 BBox3;
struct BBox3
{
	double3 min;
	double3 max;
};

BBox3 bbox_create(double3 min, double3 max);
double3 bbox_corner(BBox3 b, int corner);
BBox3 bbox_p_union(BBox3 b, double3 p);
BBox3 bbox_b_union(BBox3 b1, BBox3 b2);
BBox3 bbox_b_intersection(BBox3 b1, BBox3 b2);
bool bbox_b_overlap(BBox3 b1, BBox3 b2);
bool bbox_p_inside(double3 p, BBox3 b);
bool bbox_p_insideExclusive(double3 p, BBox3 b);
BBox3 bbox_p_expand(BBox3 b, double delta);
double3 bbox_diagonal(BBox3 b);
double bbox_surface_area(BBox3 b);
double bbox_surface_volume(BBox3 b);
int bbox_maximum_extent(BBox3 b);
double3 bbox_lerp(BBox3 b, double3 p);
double3 bbox_p_offset(BBox3 b, double3 p);

#endif
