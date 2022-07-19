#ifndef IN_OPENCL
#define IN_OPENCL
#include "bvhtree.h"
#include <opencl-c-base.h>
#endif

#include "bvhtree.h"
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

	return tmax > tmin && (tmin > 0 ? tmin : tmax) >= 0 && tmin < r.max_t;;
}

bool interactionRayTriangle(Triangle t, __constant double3 *vertices, Ray r, SurfaceInteraction *i)
{
	double dotNormalDirection, time;
	double3 C, P;
	double3 points[3] = {
		vertices[t.vertices[0]],
		vertices[t.vertices[1]],
		vertices[t.vertices[2]]
	};

	dotNormalDirection = double3_dot(t.surfaceNormal, r.direction);
	if(fabs(dotNormalDirection) < RAY_EPSILON)
		return false;

	time = -(t.det + double3_dot(t.surfaceNormal, r.origin)) / dotNormalDirection;

	if(time < r.min_t || time > r.max_t)
		return false;
	
	P = r.origin + double3_smul(time, r.direction);

	// First edge
	C = double3_cross(points[1] - points[0], P - points[0]);
	if(double3_dot(t.surfaceNormal, C) < -RAY_EPSILON)
		return false;

	// Second edge
	C = double3_cross(points[2] - points[1], P - points[1]);
	if(double3_dot(t.surfaceNormal, C) < -RAY_EPSILON)
		return false;

	// Third edge
	C = double3_cross(points[0] - points[2], P - points[2]);
	if(double3_dot(t.surfaceNormal, C) < -RAY_EPSILON)
		return false;

	i->time = time;
	i->n = t.surfaceNormal;
	i->p = P;
	
	return true;
}

bool traverseScene(Ray r, __constant double3 *vertices, __constant double3 *normals,
	__constant double2 *uv, __constant Triangle *tris,
	__constant ObjectMetadata *metadata,
	__constant BVHNode *bvh_nodes,
	SurfaceInteraction *interaction)
{
	int stack[MAX_TREE_DEPTH];
	int seen[MAX_TREE_DEPTH];
	bool interact;
	int last_node = 0;

	stack[0] = metadata->tree.root;
	seen[0] = 0;
	interact = false;
	while(last_node >= 0) {	
		if(seen[last_node] == 2 || stack[last_node] == -1)
			last_node --;
		else if(seen[last_node] == 1) {
			seen[last_node] = 2;
			last_node ++;
			seen[last_node] = 0;
			stack[last_node] = bvh_nodes[stack[last_node-1]].sons.y;
		} else if(!interactionRayBBox(bvh_nodes[stack[last_node]].bounds, r))
			last_node --;
		else if(bvh_nodes[stack[last_node]].sons.x == -1 && bvh_nodes[stack[last_node]].sons.y == -1) {
			for(int i = bvh_nodes[stack[last_node]].triangle_start; i < bvh_nodes[stack[last_node]].triangle_end; i ++) {
				if(interactionRayTriangle(tris[i], vertices, r, interaction)) {
					r.max_t = interaction->time;
					interact = true;
				}
			}
			last_node --;
		} else {
			seen[last_node] = 1;
			last_node ++;
			seen[last_node] = 0;
			stack[last_node] = bvh_nodes[stack[last_node-1]].sons.x;
		}
	}

	return interact;
}

kernel void compute_ray(
	__constant Camera *camera, __constant Material *materials,
	__constant double3 *vertices, __constant double3 *normals,
	__constant double2 *uv, __constant Triangle *tris,
	__constant ObjectMetadata *metadata,
	__constant BVHNode *bvh_nodes,
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
		.max_t = camera->max_t
	};
	
	if(!traverseScene(ray, vertices, normals, uv, tris, metadata, bvh_nodes, &interaction)) {
		write_imagei(output, on_screen_pos, default_col);
		return;
	}

	const double angle = fabs(double3_dot(interaction.n, ray.direction)); 
	const int4 col = {angle * 255, angle * 255, angle * 255, angle * 255};
	write_imagei(output, on_screen_pos, col);
} 	
