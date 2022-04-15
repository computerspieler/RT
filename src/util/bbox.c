#include "bbox.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

BBox3 bbox_create(Point min, Point max)
{
	return (BBox3) {
		.min = (Point) vec3_min(min, max),
		.max = (Point) vec3_max(min, max)
	};
}

Point bbox_corner(BBox3 b, int corner)
{
	return (Point) {
		.x = (corner & 1) ? b.min.x : b.max.x,
		.y = (corner & 2) ? b.min.x : b.max.x,
		.z = (corner & 4) ? b.min.x : b.max.x,
	};
}

BBox3 bbox_p_union(BBox3 b, Point p)
{
	return bbox_create(
		(Point) vec3_min(b.min, p),
		(Point) vec3_max(b.max, p)
	);
}

BBox3 bbox_b_union(BBox3 b1, BBox3 b2)
{
	return bbox_create(
		(Point) vec3_min(b1.min, b2.min),
		(Point) vec3_max(b1.max, b2.max)
	);
}

BBox3 bbox_b_intersection(BBox3 b1, BBox3 b2)
{
	return bbox_create(
		(Point) vec3_max(b1.min, b2.min),
		(Point) vec3_min(b1.max, b2.max)
	);
}

bool bbox_b_overlap(BBox3 b1, BBox3 b2)
{
	bool x = (b1.max.x >= b2.min.x) && (b1.min.x <= b2.max.x);
	bool y = (b1.max.y >= b2.min.y) && (b1.min.y <= b2.max.y);
	bool z = (b1.max.z >= b2.min.z) && (b1.min.z <= b2.max.z);

	return (x && y && z);
}

bool bbox_p_inside(Point p, BBox3 b)
{
	return (p.x >= b.min.x && p.x <= b.max.x &&
            p.y >= b.min.y && p.y <= b.max.y &&
            p.z >= b.min.z && p.z <= b.max.z);
}

bool bbox_p_insideExclusive(Point p, BBox3 b)
{

}