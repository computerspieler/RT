#define MAX_DEPTH_BEFORE_ROULETTE	4
#define MAX_DEPTH					20
#define RAY_EPSILON					1e-4

#include "bvhtree.h"
#include "vec3.h"
#include "bbox.h"
#include "camera.h"
#include "typedef.h"
#include "scene.h"
#include "ray.h"
#include "transform.h"
#include "scene.h"
#include "typedef.h"
#include "host.h"

#ifndef IN_OPENCL
#define IN_OPENCL
// Cette partie ne sera jamais compilé
// C'est juste pour que l'auto-complétion
// fonctionne
#define MAX_TREE_DEPTH 0
#include <opencl-c-base.h>
#include <math.h>

#endif

typedef struct Context Context;
struct Context
{
	SurfaceInteraction interaction;
	__constant Material *materials;
	__constant vec3 *vertices;
	__constant vec2 *uv;
	__constant Triangle *tris;
	__constant BVHNode *bvh_nodes;
	__constant Light *lights;
	__constant HostContext *host_ctx;

	__global int *selected_triangle_ptr;

	int selected_triangle;
	bool selectedByCursor;
	uint state;
	int2 on_screen_pos;
};

typedef struct PathElement PathElement;
struct PathElement
{
	vec3 p;
	vec3 wi, wo;
	Float weight;
    int triangle_id;
	bool normalReversed;
	Transform worldObject;
	Float current_IOR;
};

uint randInt(Context *ctx)
{
	uint x = ctx->state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return ctx->state = x;
}

Float randFloat(Context *ctx)
{
	return fabs((Float) (randInt(ctx)) / (Float) (0xFFFFFFFF));
}

bool material_is_specular(Material m)
{
	switch(m.type) {
	case MATERIAL_NONE: return false;
    case MATERIAL_LAMBERTIAN: return false;
    case MATERIAL_MIRROR: return true;
    case MATERIAL_GLASS: return true;
	}
}

vec3 UniformSphereSampling(Context *ctx)
{
	const Float u1 = 2. * randFloat(ctx) - 1.;
	const Float u2 = randFloat(ctx);

    const Float phi = 2. * M_PI * u2;
    const Float r = sqrt(1. - u1 * u1);
	vec3 output;

	output.x = cos(phi) * r;
	output.y = sin(phi) * r;
	output.z = u1;

    return output;
}

// Source: https://www.rorydriscoll.com/2009/01/07/better-sampling/
vec3 UniformSampleHemisphere(Context *ctx)
{
	const Float u1 = randFloat(ctx);
	const Float u2 = randFloat(ctx);

    const Float phi = 2. * M_PI * u2;
    const Float r = sqrt(1. - u1 * u1);
	vec3 output;

	output.x = cos(phi) * r;
	output.y = sin(phi) * r;
	output.z = u1;

    return output;
}

// Source: https://www.rorydriscoll.com/2009/01/07/better-sampling/
vec3 CosineSampleHemisphere(Context *ctx)
{
	const Float u1 = randFloat(ctx);
	const Float u2 = randFloat(ctx);

    const Float theta = 2. * M_PI * u2;
    const Float r = sqrt(u1);
	vec3 output;

	output.x = cos(theta) * r;
	output.y = sin(theta) * r;
	output.z = sqrt(fmax(0., 1. - u1));

    return output;
}

bool interactionRayBBox(BBox3 b, Ray r, Float *time)
{
	vec3 t1, t2;
	Float tmin, tmax;

	t1 = (b.min - r.origin) / r.direction;
	t2 = (b.max - r.origin) / r.direction;

	tmin = max(
		max(
			(Float) min(t1.x, t2.x),
			(Float) min(t1.y, t2.y)
		),  (Float) min(t1.z, t2.z)
	);

	tmax = min(
		min(
			(Float) max(t1.x, t2.x),
			(Float) max(t1.y, t2.y)
		),  (Float) max(t1.z, t2.z)
	);

	if(time) *time = tmin > 0 ? tmin : tmax;
	return tmax >= tmin && (tmin > 0 ? tmin : tmax) >= 0 && tmin <= r.max_t;
}

bool interactionRaySphere(Ray r, vec3 center, Float radius, Float *time)
{
	vec3 L = r.origin - center;
	Float a = vec3_norm_2(r.direction);
	Float b = vec3_dot(r.direction, L);
	Float c = vec3_norm_2(L) - radius*radius;

	Float discr = b*b - a*c;
	if(discr <= 0) return false;

	discr = sqrt(discr);
	Float t0 = (- b - discr) / a;
	Float t1 = (- b + discr) / a;

	if(t0 > t1) {
		Float tmp = t0;
		t0 = t1;
		t1 = tmp;
	}

	if(t0 < 0)
		t0 = t1;

	if(time) *time = t0;
	return t0 > r.min_t && t0 < r.max_t;
}

bool interactionRayLight(Ray r, Light l, Float *time)
{
	BBox3 b;
	vec3 tmp;
	switch(l.type) {
	case LIGHT_POINT:
		return interactionRaySphere(r, l.pos, 1e-2, time);
	
	case LIGHT_AREA:
		tmp =
			(Float) (.5)   * l.area_width   * l.area_tan   +
			(Float) (.5)   * l.area_height  * l.area_bitan +
			(Float) (RAY_EPSILON) * l.area_n;
		b = (BBox3) {
			.min = l.pos - tmp,
			.max = l.pos + tmp
		};

		return interactionRayBBox(b, r, time);
	}
}

bool interactionRayTriangle(uint t_id, Ray r, Context *ctx)
{
	vec3 d;
	vec3 e1, e2;
	vec3 s1, s2;
	vec3 p1;
	vec2 uv0, uv1, uv2;
	Float div, invDiv;
	Float b0, b1, b2;
	Float t;

	p1 = ctx->vertices[ctx->tris[t_id].vertices[0]];

	e1 = ctx->vertices[ctx->tris[t_id].vertices[1]] - p1;
	e2 = ctx->vertices[ctx->tris[t_id].vertices[2]] - p1;
	s1 = vec3_cross(r.direction, e2);
	div = vec3_dot(s1, e1);
	if(div == 0.)
		return false;
	invDiv = 1. / div;

	d = r.origin - p1;
	b1 = vec3_dot(d, s1) * invDiv;
	if(b1 < 0. || b1 > 1.)
		return false;

	s2 = vec3_cross(d, e1);
	b2 = vec3_dot(r.direction, s2) * invDiv;
	if(b2 < 0. || b1 + b2 > 1.)
		return false;

	b0 = 1. - b1 - b2;

	t = vec3_dot(e2, s2) * invDiv;
	if(t < r.min_t || t > r.max_t)
		return false;

	uv0 = ctx->uv[ctx->tris->uv[0]];
	uv1 = ctx->uv[ctx->tris->uv[1]];
	uv2 = ctx->uv[ctx->tris->uv[2]];

	ctx->interaction = (SurfaceInteraction) {
		.m = ctx->materials[ctx->tris[t_id].material],
		.time = t,
		.p = r.origin + t * r.direction,
		
		.n = ctx->tris[t_id].normal,
		.dpdu = ctx->tris[t_id].dpdu,
		.dpdv = ctx->tris[t_id].dpdv,

		.uv = b0 * uv0 + b1 * uv1 + b2 * uv2
	};

	return true;
}

bool traverseScene(Ray r, Context *ctx)
{
	int i;
	bool interact;
	interact = false;
#if 0
	for(i = 0; i < ctx->host_ctx->scene_meta.triangles_count; i ++) {
		if(interactionRayTriangle(i, r, ctx)) {
			r.max_t = ctx->interaction.time;
			ctx->interaction.triangle_id = i;
			interact = true;
		}
	}
#else
	int stack[MAX_TREE_DEPTH + 1];
	int seen [MAX_TREE_DEPTH + 1];
	int last_node = 0;

	seen[0] = 0;
	stack[0] = ctx->host_ctx->scene_meta.tree.root;

	while(last_node >= 0) {	
		if(seen[last_node] == 2 || stack[last_node] == -1)
			last_node --;
		else if(seen[last_node] == 1) {
			seen[last_node] = 2;
			last_node ++;
			seen[last_node] = 0;
			stack[last_node] = ctx->bvh_nodes[stack[last_node-1]].sons.y;
		}
		else if(!interactionRayBBox(ctx->bvh_nodes[stack[last_node]].bounds, r, NULL))
			last_node --;
		
		// Si c'est une feuille
		else if(ctx->bvh_nodes[stack[last_node]].sons.x == -1 &&
				ctx->bvh_nodes[stack[last_node]].sons.y == -1) {
			for(i = ctx->bvh_nodes[stack[last_node]].triangle_start;
				i < ctx->bvh_nodes[stack[last_node]].triangle_end;
				i ++) {
				if(interactionRayTriangle(i, r, ctx)) {
					r.max_t = ctx->interaction.time;
					ctx->interaction.triangle_id = i;
					interact = true;
				}
			}
			last_node --;
		}
		else {
			seen[last_node] = 1;
			last_node ++;
			seen[last_node] = 0;
			stack[last_node] = ctx->bvh_nodes[stack[last_node-1]].sons.x;
		}
	}
#endif
	return interact;
}

vec3 getPointFromLightSource(Context *ctx, int light_id) {
	vec3 output = ctx->lights[light_id].pos;

	switch(ctx->lights[light_id].type) {
	case LIGHT_POINT: break;
	case LIGHT_AREA:
		output +=
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_width ) * ctx->lights[light_id].area_tan +
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_height) * ctx->lights[light_id].area_bitan;
		break;
	}

	return output;
}

Float light_pdf(Context *ctx, int light_id) {
	switch(ctx->lights[light_id].type) {
	case LIGHT_POINT:
		return .25 * M_1_PI; // 1/(4*pi)
	case LIGHT_AREA:
		return .25 * M_1_PI / (ctx->lights[light_id].area_width * ctx->lights[light_id].area_height);
	}
}

Float light_l(Context *ctx, int light_id, vec3 ref_point, vec3 light_point, bool check_for_collision)
{
	Ray ray;
	Float output;
	vec3 dist, wo;

	if(!check_for_collision)
		return ctx->lights[light_id].I / vec3_norm_2(light_point - ref_point);

	dist = light_point - ref_point;
	wo = vec3_normalize(dist);

	ray = (Ray) {
		.origin = ref_point,
		.max_t = vec3_norm(dist),
		.min_t = RAY_EPSILON,
		.direction = wo
	};

	if(ctx->selectedByCursor)
		printf("Light to last point's distance: " FLOAT_FMT "\n", ray.max_t);
	
	output = 0.;
	if(!traverseScene(ray, ctx)) {
		if(ctx->selectedByCursor)
			printf("No collision ahead\n");
	
		output = ctx->lights[light_id].I
			/ (vec3_norm_2(dist)) // * light_pdf(ctx, light_id))
			;
	}
	else if(ctx->selectedByCursor)
		printf("Time of collision: " FLOAT_FMT "; with the triangle %d\n",
			ctx->interaction.time, ctx->interaction.triangle_id);
	
	return output;
}

Float material_compute_brdf(Material *m, vec3 wi, vec3 wo, int lambda, Float current_IOR, Float *new_IOR)
{
    vec3 valid_wo;
	Float IOR;
	Float sin_theta_refracted, theta_refracted;
	
    switch(m->type) {
        case MATERIAL_NONE: return 0;
        default: return 0;

        case MATERIAL_LAMBERTIAN:
            return m->l.R * M_1_PI;
        
		case MATERIAL_MIRROR:
			if(wi.x == wo.x && wi.y == wo.y && wi.z == -wo.z)
            	return 1;
			else
				return 0;
		
        case MATERIAL_GLASS:
			return 1;
    }
}

vec3 material_sample_wo(Material *m, vec3 wi, Context *ctx, Float current_IOR, Float *pdf)
{
	vec3 wo, mirror_wo, refracted_wo;
	Float IOR, tmp;
	Float T_TE;
	Float sin_2_theta_incident;
	Float sin_2_theta_refracted, sin_theta_refracted;

	wo = VEC3(0, 0, 0);

    switch(m->type) {
        default:
        case MATERIAL_NONE:
			if(pdf) *pdf = 1;
			return wo;
		
        case MATERIAL_GLASS:
			IOR = current_IOR == m->g.IOR ? 1 : m->g.IOR;
			sin_2_theta_incident = wi.x * wi.x + wi.y * wi.y;

			mirror_wo = wi;
			mirror_wo.z *= -1;

			tmp = current_IOR / IOR; tmp *= tmp;
			sin_2_theta_refracted = tmp * sin_2_theta_incident;

			if(sin_2_theta_refracted > 1) {
				if(pdf) *pdf = 1;
				return mirror_wo;
			}
			
			sin_theta_refracted = sqrt(sin_2_theta_refracted);

			refracted_wo.x = wi.x * sin_theta_refracted;
			refracted_wo.y = wi.y * sin_theta_refracted;
			refracted_wo.z = wi.z;

			refracted_wo = vec3_normalize(refracted_wo);

			T_TE = 2. * current_IOR * fabs(wi.z);
			T_TE /= current_IOR * fabs(wi.z) + IOR * fabs(refracted_wo.z);
#if PRINT_DEBUG>=3
			if(ctx->selectedByCursor) {
				printf("Here we go ! " FLOAT_FMT "\n", T_TE);
			}
#endif

			if(randFloat(ctx) < T_TE) {
				if(pdf) *pdf = T_TE;
				return refracted_wo;
			} else {
				if(pdf) *pdf = 1. - T_TE;
				return mirror_wo;
			}

        case MATERIAL_MIRROR:
			wo = wi;
            wo.z *= -1;
			if(pdf) *pdf = 1;
			return wo;
		
        case MATERIAL_LAMBERTIAN:
			wo = CosineSampleHemisphere(ctx);
			wo.z *= -1;
			if(pdf) *pdf = fabs(wo.z) * M_1_PI;
			return wo;
    }
}

Float compute_simple_ray(Ray ray, Context *ctx)
{
	int i, j;
	Float x, y, z;
	vec3 n;
	//vec3 dpdu, dpdv;
	
	if(!traverseScene(ray, ctx)) {
		if(ctx->selectedByCursor) {
			ctx->selected_triangle = -1;
			*(ctx->selected_triangle_ptr) = -1;
		}

#if PRINT_DEBUG>=3
		if(ctx->selectedByCursor)
			printf("Out !\n");
#endif
		return 0;
	}

	n = ctx->interaction.n;

	if(ctx->selectedByCursor) {
		ctx->selected_triangle = ctx->interaction.triangle_id;
		*(ctx->selected_triangle_ptr) = ctx->interaction.triangle_id;
	}

#if PRINT_DEBUG>=3
	if(ctx->selectedByCursor) {
		printf("Normal: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
			ctx->interaction.n.x, ctx->interaction.n.y, ctx->interaction.n.z);
		printf("DPDU: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
			ctx->interaction.dpdu.x, ctx->interaction.dpdu.y, ctx->interaction.dpdu.z);
		printf("DPDV: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
			ctx->interaction.dpdv.x, ctx->interaction.dpdv.y, ctx->interaction.dpdv.z);
		printf("Hit: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
			ctx->interaction.p.x, ctx->interaction.p.y, ctx->interaction.p.z);
		printf("UV: " FLOAT_FMT " " FLOAT_FMT "\n",
			ctx->interaction.uv.x, ctx->interaction.uv.y);
		printf("Time : " FLOAT_FMT "\n",
			ctx->interaction.time);
	}
#endif

	if(ctx->selected_triangle == ctx->interaction.triangle_id)
		return 1.;
	else
		return fabs(vec3_dot(ray.direction, n));
}


int compute_ray(Ray ray, Context *ctx, PathElement *path)
{
	vec3 wi, wo;
	int depth;
	Transform worldObject;
	SurfaceInteraction interaction;

	bool addLastPoint;
	int i, j;
    Float weight;
	Float x, y, z;
	Float rrFactor;
	Float rrStopProbability;
	Float current_IOR, new_IOR;
	Float pdf;

    weight = 1;
	current_IOR = 1;

	path[0] = (PathElement){
		.p = ray.origin,
		.weight = 1,
		.wi = VEC3(0, 0, 0),
		.wo = VEC3(0, 0, 0),
		.normalReversed = false,
        .triangle_id = -1,
		.current_IOR = current_IOR
	};

	for(depth = 1; weight > 0 && depth <= MAX_DEPTH; depth ++) {
		rrFactor = 1.0;
		if(depth > MAX_DEPTH_BEFORE_ROULETTE) {
			rrStopProbability = 0.9;
			Float d = randFloat(ctx);
            if(d <= rrStopProbability) {
#if PRINT_DEBUG>=3
				if(ctx->selectedByCursor)
					printf("Roulette decided: " FLOAT_FMT "\n", d);
#endif
				break;
			}
			rrFactor = 1.0 / (1.0 - rrStopProbability);
		}

#if PRINT_DEBUG>=3
		if(ctx->selectedByCursor) {
			printf("Depth: %d\n", depth);
			printf("Wi: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				ray.direction.x, ray.direction.y, ray.direction.z);
			printf("Origin: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				ray.origin.x, ray.origin.y, ray.origin.z);
		}
#endif

		if(!traverseScene(ray, ctx)) {
#if PRINT_DEBUG>=3
			if(ctx->selectedByCursor)
				printf("Out !\n");
#endif
			break;
		}

		interaction = ctx->interaction;
		if(vec3_dot(ray.direction, ctx->interaction.n) < 0) {
			interaction.n *= -1;
			interaction.dpdu *= -1;
		}

		// Prepare the matrix
		worldObject.mat.v4[0][0] = interaction.dpdu.x;
		worldObject.mat.v4[0][1] = interaction.dpdv.x;
		worldObject.mat.v4[0][2] = interaction.n.x;
		worldObject.mat.v4[0][3] = 0;

		worldObject.mat.v4[1][0] = interaction.dpdu.y;
		worldObject.mat.v4[1][1] = interaction.dpdv.y;
		worldObject.mat.v4[1][2] = interaction.n.y;
		worldObject.mat.v4[1][3] = 0;

		worldObject.mat.v4[2][0] = interaction.dpdu.z;
		worldObject.mat.v4[2][1] = interaction.dpdv.z;
		worldObject.mat.v4[2][2] = interaction.n.z;
		worldObject.mat.v4[2][3] = 0;

		worldObject.mat.v4[3][0] = 0;
		worldObject.mat.v4[3][1] = 0;
		worldObject.mat.v4[3][2] = 0;
		worldObject.mat.v4[3][3] = 1;

		/*
			On profite du fait que la base (dpdu, dpdv, n) soit une
			base orthonormée, pour optimiser le code
		*/
		//matrix_inverse(&worldObject.matInv, worldObject.mat);
		matrix_transpose(&worldObject.matInv, worldObject.mat);

		weight *= rrFactor;

		// Compute the ray's next direction
		wi = vec3_normalize(transform_apply_vector(ray.direction, transform_inverse(worldObject)));
		wo = material_sample_wo(&interaction.m, wi, ctx, current_IOR, &pdf);
		
		if(wo.x == 0 && wo.y == 0 && wo.z == 0) {
#if PRINT_DEBUG>=3
			if(ctx->selectedByCursor)
				printf("No more WO !\n");
#endif
			break;
		}
		
		weight = weight * fabs(wo.z)
			* material_compute_brdf(&interaction.m, wi, wo, ctx->host_ctx->lambda, current_IOR, &new_IOR)
			/ pdf;
	
#if PRINT_DEBUG>=3
		if(ctx->selectedByCursor) {
			printf("Normal: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				interaction.n.x, interaction.n.y, interaction.n.z);
			printf("DPDU: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				interaction.dpdu.x, interaction.dpdu.y, interaction.dpdu.z);
			printf("DPDV: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				interaction.dpdv.x, interaction.dpdv.y, interaction.dpdv.z);
			printf("Hit: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				interaction.p.x, interaction.p.y, interaction.p.z);
			printf("Wi: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				wi.x, wi.y, wi.z);
			printf("Wo: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
				wo.x, wo.y, wo.z);
			printf("UV: " FLOAT_FMT " " FLOAT_FMT "\n",
				interaction.uv.x, interaction.uv.y);
			printf("BRDF: " FLOAT_FMT "\n", material_compute_brdf(&interaction.m, wi, wo, ctx->host_ctx->lambda, current_IOR, NULL));
			printf("Time : " FLOAT_FMT "\nWeight : " FLOAT_FMT "\nTriangle : %d\n",
				interaction.time, weight, interaction.triangle_id);
			matrix_print(worldObject.mat);

		}
#endif

		path[depth] = (PathElement){
			.p = interaction.p,
			.weight = weight,
			.wi = wi,
			.wo = wo,
			.normalReversed = vec3_dot(ray.direction, ctx->interaction.n) < 0,
            .triangle_id = interaction.triangle_id,
			.worldObject = worldObject,
			.current_IOR = current_IOR
		};

		// Calcul du rayon de la prochaine étape
		ray.origin = path[depth].p;
		ray.max_t = ctx->host_ctx->camera.max_t;
		ray.direction = vec3_normalize(transform_apply_vector(wo, worldObject));
		current_IOR = new_IOR;
	}

	return depth - 1;
}

kernel void compute_pixel(
	__constant Material *materials,
	__constant vec3 *vertices,
	__constant vec2 *uv,
	__constant Triangle *tris,
	__constant BVHNode *bvh_nodes,
	__constant HostContext *host_ctx,
	__global int *selected_triangle,
	__constant Light *lights,
	__global Float *raw,
	__write_only image2d_t output,
	__global Float *surface_irradiance
)
{
	int i, j;
	const int2 on_screen_pos = {get_global_id(0), get_global_id(1)};

	Context ctx = {
		.interaction = (SurfaceInteraction) {0},
		.materials = materials,
		.vertices = vertices,
		.uv = uv,
		.tris = tris,
		.bvh_nodes = bvh_nodes,
		.host_ctx = host_ctx,
		.lights = lights,
		.on_screen_pos = on_screen_pos,
		.selected_triangle = *selected_triangle,
		.selected_triangle_ptr = selected_triangle,
		.selectedByCursor =
			on_screen_pos.x == host_ctx->mouse_pos.x &&
			on_screen_pos.y == host_ctx->mouse_pos.y,
		.state =
			host_ctx->seed + on_screen_pos.x * on_screen_pos.y
	};

	if(ctx.state % 2 == 0)
		ctx.state ++;

	vec2 normalized_screen_pos = (vec2) (
		((Float) (2 * on_screen_pos.x) / (Float) host_ctx->camera.viewport.x) - 1.0f,
		((Float) (2 * on_screen_pos.y) / (Float) host_ctx->camera.viewport.y) - 1.0f
	);
	Float l = tan(host_ctx->camera.fov) * host_ctx->camera.near;
	vec3 tmp = (vec3)(
		normalized_screen_pos.x * l,
		normalized_screen_pos.y * l,
		host_ctx->camera.near
	);
	tmp = transform_apply_vector(tmp, host_ctx->camera.rotation_transform);
	Ray ray = (Ray){
		.origin = tmp + host_ctx->camera.pos,
		.direction = vec3_normalize(tmp),
		.min_t = 1e-5,
		.max_t = host_ctx->camera.max_t
	};
	
	const int index = on_screen_pos.y + host_ctx->camera.viewport.y * on_screen_pos.x;
	Float grayscale;

	if(ctx.host_ctx->simpleView) {
		for(i = 0; i < host_ctx->scene_meta.lights_count; i ++)
			if(interactionRayLight(ray, ctx.lights[i], NULL))
				break;
		if(i != host_ctx->scene_meta.lights_count) {
			// On colorie la face en noir si elle ne nous fait pas face
			if(ctx.lights[i].type == LIGHT_AREA && vec3_dot(ray.direction, ctx.lights[i].area_n) < 0)
				grayscale = 0;
			else
				grayscale = 255;
		} else
			grayscale = 127 * compute_simple_ray(ray, &ctx);
		write_imagei(output, on_screen_pos, INT4(grayscale, grayscale, grayscale, 0));
	} else {
		Float tmp, totalVal, weight, light_l_pow;
		int cam_path_length = 0;
		PathElement cam_path[MAX_DEPTH + 1];
		
#if PRINT_DEBUG>=2
		if(ctx.selectedByCursor)
			printf("Computing camera's path\n");
#endif

		cam_path_length = compute_ray(ray, &ctx, cam_path);

		for(i = 1; i <= cam_path_length; i ++) {	
			int light_id = randFloat(&ctx) * host_ctx->scene_meta.lights_count;		
			if(ctx.selectedByCursor)
				printf("Light %d for %d\n", light_id, i);
			
			Material mat = materials[tris[cam_path[i].triangle_id].material];
			weight = 0;

			if(!material_is_specular(mat)) {
				vec3 light_point = getPointFromLightSource(&ctx, light_id);
				light_l_pow = light_l(&ctx, light_id, cam_path[i].p, light_point, true);
				
				if(light_l_pow == 0)
					continue;

				vec3 wo = transform_apply_vector(
					vec3_normalize(light_point - cam_path[i].p),
					transform_inverse(cam_path[i].worldObject)
				);
				weight = cam_path[i-1].weight
					* material_compute_brdf(&mat, cam_path[i].wi, wo, host_ctx->lambda, cam_path[i].current_IOR, NULL)
					* fabs(wo.z);
				
				if(ctx.selectedByCursor) {
					printf("Cam's Point: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
						cam_path[i].p.x, cam_path[i].p.y, cam_path[i].p.z);
					printf("Light's Point: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
						light_point.x, light_point.y, light_point.z);
					printf("Light's power : " FLOAT_FMT "; Path's weight: " FLOAT_FMT "\n", light_l_pow, weight);
					printf("Has the normal been reversed ? %s\n", cam_path[i].normalReversed ? "Yes" : "No");
				}
			} else {
				Ray r = (Ray) {
					.origin = cam_path[i].p,
					.max_t =
						i < cam_path_length - 1 ?
							vec3_norm(cam_path[i+1].p - cam_path[i].p) :
							host_ctx->camera.max_t,
					.min_t = RAY_EPSILON,
					.direction = transform_apply_vector(
						material_sample_wo(&mat, cam_path[i].wi, &ctx, cam_path[i].current_IOR, NULL),
						cam_path[i].worldObject
					)
				};

				if(ctx.selectedByCursor) {
					printf("Cam's Point: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n",
						cam_path[i].p.x, cam_path[i].p.y, cam_path[i].p.z);
				}

				Float time;
				if(interactionRayLight(r, ctx.lights[light_id], &time)) {
					weight = cam_path[i].weight;
					light_l_pow = light_l(&ctx, light_id, cam_path[i].p, cam_path[i].p + time * r.direction, false);

					if(ctx.selectedByCursor) {
						printf("Went through : " FLOAT_FMT "; " FLOAT_FMT " !\n", time, light_l_pow);
					}
				}
			}

			if(ctx.selectedByCursor) {
				printf("Added : " FLOAT_FMT "; " FLOAT_FMT " !\n", weight, light_l_pow);
			}
			totalVal += weight * light_l_pow;
		}
		if(ctx.selectedByCursor) {
			printf("Total : " FLOAT_FMT " !\n", totalVal);
		}

		if(host_ctx->sampleCount == 1)
			raw[index] = 0;
		if(cam_path_length > 0)
			raw[index] += totalVal / cam_path_length;

#if PRINT_DEBUG>=1
		if(ctx.selectedByCursor)
			printf("Mean weight : " FLOAT_FMT "; Camera's path's length: %d; Selected triangle: %d\n",
				raw[index] / host_ctx->sampleCount, cam_path_length, ctx.selected_triangle);
#endif

		tmp = raw[index] / host_ctx->sampleCount;
        if(tmp < 0)
			write_imagei(output, on_screen_pos, INT4(0, 255, 0, 0));
		else {
			grayscale = 255 * log(tmp + 1) / host_ctx->max_threshold;
            if(grayscale >= 255) grayscale = 255;
			write_imagei(output, on_screen_pos, INT4(grayscale, grayscale, grayscale, 0));
		} 
	}
}

