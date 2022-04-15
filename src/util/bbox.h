#ifndef _BBOX_H_
#define _BBOX_H_

#include <math.h>
#include <stdbool.h>

#include "point.h"

typedef struct BBox3 BBox3;
struct BBox3
{
	Point min;
	Point max;
};

const Point default_min = (Point) {
	.x = INFINITY,
	.y = INFINITY,
	.z = INFINITY
};

const Point default_max = (Point) {
	.x = -INFINITY,
	.y = -INFINITY,
	.z = -INFINITY
};

BBox3 bbox_create(Point min, Point max);
Point bbox_corner(BBox3, int corner);
BBox3 bbox_p_union(BBox3, Point);
BBox3 bbox_b_union(BBox3, BBox3);

#endif
