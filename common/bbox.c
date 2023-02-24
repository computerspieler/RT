#include "typedef.h"

#include "vec3.h"
#include "scene.h"
#include "bbox.h"

BBox3 bbox_create(vec3 min, vec3 max)
{
	return (BBox3) {
		.min = vec3_min(min, max),
		.max = vec3_max(min, max)
	};
}

vec3 bbox_corner(BBox3 b, int corner)
{
	vec3 f;
	f.x = (corner & 1) ? b.min.x : b.max.x;
	f.y = (corner & 2) ? b.min.y : b.max.y;
	f.z = (corner & 4) ? b.min.z : b.max.z;
	return f;
}

vec3 bbox_center(BBox3 b)
{
	return vec3_smul(0.5, vec3_add(b.max, b.min));
}

BBox3 bbox_p_union(BBox3 b, vec3 p)
{
	return bbox_create(
		(vec3) vec3_min(b.min, p),
		(vec3) vec3_max(b.max, p)
	);
}

BBox3 bbox_b_union(BBox3 b1, BBox3 b2)
{
	return bbox_create(
		(vec3) vec3_min(b1.min, b2.min),
		(vec3) vec3_max(b1.max, b2.max)
	);
}

BBox3 bbox_b_intersection(BBox3 b1, BBox3 b2)
{
	return bbox_create(
		(vec3) vec3_max(b1.min, b2.min),
		(vec3) vec3_min(b1.max, b2.max)
	);
}

bool bbox_b_overlap(BBox3 b1, BBox3 b2)
{
	bool x = (b1.max.x >= b2.min.x) && (b1.min.x <= b2.max.x);
	bool y = (b1.max.y >= b2.min.y) && (b1.min.y <= b2.max.y);
	bool z = (b1.max.z >= b2.min.z) && (b1.min.z <= b2.max.z);

	return (x && y && z);
}

bool bbox_p_inside(vec3 p, BBox3 b)
{
	return (p.x >= b.min.x && p.x <= b.max.x &&
			p.y >= b.min.y && p.y <= b.max.y &&
			p.z >= b.min.z && p.z <= b.max.z);
}

bool bbox_p_insideExclusive(vec3 p, BBox3 b)
{
	return (p.x >= b.min.x && p.x < b.max.x &&
            p.y >= b.min.y && p.y < b.max.y &&
            p.z >= b.min.z && p.z < b.max.z);
}

BBox3 bbox_p_expand(BBox3 b, Float delta)
{
	vec3 v_delta;
	v_delta.x = v_delta.y = v_delta.z = delta;

	return bbox_create(
		vec3_diff(b.min, v_delta),
		vec3_add (b.max, v_delta)
	);
}

vec3 bbox_diagonal(BBox3 b)
{
	return vec3_diff(b.max, b.min);
}

Float bbox_surface_area(BBox3 b)
{
	vec3 diag = bbox_diagonal(b);
	return 2 * vec3_norm_2(diag);
}

Float bbox_surface_volume(BBox3 b)
{
	vec3 diag = bbox_diagonal(b);
	return diag.x * diag.y * diag.z;
}

int bbox_maximum_extent(BBox3 b)
{
	vec3 diag = bbox_diagonal(b);
	if(diag.x > diag.y && diag.x > diag.z)
		return 0;
	else if(diag.y > diag.z)
		return 1;
	else
		return 2;
}

vec3 bbox_lerp(BBox3 b, vec3 p)
{
	return vec3_lerp(b.min, b.max, p);
}

vec3 bbox_p_offset(BBox3 b, vec3 p)
{
	vec3 o = vec3_diff(p, b.min);

	if(b.max.x > b.min.x) o.x /= b.max.x - b.min.x;
	if(b.max.y > b.min.y) o.y /= b.max.y - b.min.y;
	if(b.max.z > b.min.z) o.z /= b.max.z - b.min.z;

	return o;
}
