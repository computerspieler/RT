#ifndef IN_OPENCL
#define IN_OPENCL
#include <opencl-c-base.h>
#endif

#include "double3.h"
#include "bbox.h"
#include "camera.h"
#include "typedef.h"
#include "scene.h"
#include "ray.h"
#include "transform.h"
#include "object.h"

bool interactionRayBBox(BBox3 b, Ray r)
{
	double3 t1, t2;
	double tmin, tmax;

	t1 = (b.min - r.origin) / r.direction;
	t2 = (b.max - r.origin) / r.direction;

	tmin = max(
		max(
			(double) min(t1.x, t2.x),
			(double) min(t1.y, t2.y)
		), (double) min(t1.z, t2.z)
	);

	tmax = min(
		min(
			(double) max(t1.x, t2.x),
			(double) max(t1.y, t2.y)
		), (double) max(t1.z, t2.z)
	);

	return tmax > tmin && (tmin > 0 ? tmin : tmax) >= 0 && tmin < r.max_t;
}

bool interactionRayTriangle(Triangle t, __global double3 *vertices, Ray r, SurfaceInteraction *i)
{
	double dotNormalDirection, time, det;
	double3 C, P, surfaceNormal;
	double3 points[3] = {
		vertices[t.vertices[0]],
		vertices[t.vertices[1]],
		vertices[t.vertices[2]]
	};
	
	surfaceNormal = double3_cross(double3_diff(points[1], points[0]), double3_diff(points[2], points[0]));
    surfaceNormal = double3_normalize(surfaceNormal);
	det = -double3_dot(surfaceNormal, points[0]);

	dotNormalDirection = double3_dot(surfaceNormal, r.direction);
	if(fabs(dotNormalDirection) < RAY_EPSILON)
		return false;

	time = -(det + double3_dot(surfaceNormal, r.origin)) / dotNormalDirection;

	if(time < r.min_t || time > r.max_t)
		return false;
	
	P = r.origin + double3_smul(time, r.direction);

	// First edge
	C = double3_cross(points[1] - points[0], P - points[0]);
	if(double3_dot(surfaceNormal, C) < -RAY_EPSILON)
		return false;

	// Second edge
	C = double3_cross(points[2] - points[1], P - points[1]);
	if(double3_dot(surfaceNormal, C) < -RAY_EPSILON)
		return false;

	// Third edge
	C = double3_cross(points[0] - points[2], P - points[2]);
	if(double3_dot(surfaceNormal, C) < -RAY_EPSILON)
		return false;

	i->time = time;
	i->n = surfaceNormal;
	i->p = P;
	
	return true;
}

kernel void compute_ray(
	__global Camera *camera, __global Material *materials,
	__global double3 *vertices, __global double3 *normals,
	__global double3 *uv, __global Triangle *tris,
	__global ObjectMetadata *metadata,
	__write_only image2d_t output
)
{
	// BGR
	const int4 default_col = {239, 149, 102, 0};
	SurfaceInteraction interaction, nearestInteraction;
	int2 on_screen_pos = {get_global_id(0), get_global_id(1)};

	double2 normalized_screen_pos = (double2) (
		((double) (2 * on_screen_pos.x) / (double) camera->viewport.x) - 1.0f,
		((double) (2 * on_screen_pos.y) / (double) camera->viewport.y) - 1.0f
	);
	double l = tan(camera->fov) * camera->near;
	double3 tmp = (double3)(
		normalized_screen_pos.x * l,
		normalized_screen_pos.y * l,
		camera->near
	);
	Ray ray = (Ray){
		.origin = tmp + camera->pos,
		.direction = double3_normalize(tmp),
		.min_t = 0,
		.max_t = 10
	};

	bool inside;
	
	inside = interactionRayBBox(metadata->bounds, ray);

	if(!inside) {
		write_imagei(output, on_screen_pos, default_col);
		return;
	}

	int i;
	nearestInteraction.time = -1;
	for(i = 0; i < metadata->triangles_count; i++) {
		if(interactionRayTriangle(tris[i], vertices, ray, &interaction)) {
			nearestInteraction = interaction;
			ray.max_t = nearestInteraction.time;
		}
	}

	if(nearestInteraction.time < 0) {
		write_imagei(output, on_screen_pos, default_col);
		return;
	}
	
	const double angle = fabs(double3_dot(nearestInteraction.n, ray.direction)); 
	const int4 col = {angle * 255, angle * 255, angle * 255, angle * 255};
	write_imagei(output, on_screen_pos, col);
} 	
