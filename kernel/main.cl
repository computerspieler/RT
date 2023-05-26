#include "bvhtree.h"
#include "vec3.h"
#include "bbox.h"
#include "typedef.h"
#include "scene.h"
#include "ray.h"
#include "transform.h"
#include "scene.h"
#include "typedef.h"
#include "host.h"
#include "sample.h"

#ifndef IN_OPENCL
#define IN_OPENCL
// Cette partie ne sera jamais compilé
// C'est juste pour que l'auto-complétion
// fonctionne
#include <opencl-c-base.h>
#include <math.h>

#endif

#include "tinymt32.h"

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

	bool verbose;
	tinymt32_t twister;
};

typedef struct PathElement PathElement;
struct PathElement
{
	vec3 p;
	vec3 wi, wo;
	Float weight;
	Float pdf;
    int triangle_id;
	bool normalReversed;
	Transform worldObject;
	Float current_IOR;
};

uint randInt(Context *ctx)
{
	return tinymt32_generate_uint32(&ctx->twister);
}

Float randFloat(Context *ctx)
{
	return tinymt32_generate_32double(&ctx->twister);
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

bool interactionRayTriangle(uint t_id, Ray r, Context *ctx)
{
	vec3 d;
	vec3 e1, e2;
	vec3 s1, s2;
	vec3 p1;
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

	ctx->interaction = (SurfaceInteraction) {
		.m = ctx->materials[ctx->tris[t_id].material],
		.time = t,
		.p = r.origin + t * r.direction,
		
		.n = ctx->tris[t_id].normal,
		.dpdu = ctx->tris[t_id].dpdu,
		.dpdv = ctx->tris[t_id].dpdv,

		.uv = b0 * ctx->uv[ctx->tris->uv[0]]
			+ b1 * ctx->uv[ctx->tris->uv[1]]
			+ b2 * ctx->uv[ctx->tris->uv[2]]
	};

	return true;
}

bool traverseScene(Ray r, Context *ctx)
{
	int i;
	bool interact;
	interact = false;

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

	return interact;
}

Float material_compute_brdf(Material *m, vec3 wi, vec3 wo, Float current_IOR, Float *new_IOR)
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
			if(pdf) *pdf = wo.z * M_1_PI;
			wo.z *= -1;
			return wo;
    }
}

int calcul_chemin(Context *ctx, Ray ray, PathElement *chemin, bool *went_out_of_bounds)
{
	vec3 wi, wo;
	int depth;
	Transform worldObject;
	SurfaceInteraction interaction;

	bool addLastPoint;
	int i, j;
	Float x, y, z;
	Float rrFactor;
	Float rrStopProbability;
	Float current_IOR, new_IOR;
	Float pdf;

	current_IOR = 1;

	chemin[0] = (PathElement){
		.pdf = 1,
		.p = ray.origin,
		.weight = 1,
		.wi = VEC3(0, 0, 0),
		.wo = VEC3(0, 0, 0),
		.normalReversed = false,
        .triangle_id = -1,
		.current_IOR = current_IOR
	};

	// On prépare les valeurs qui ne changeront pas
	worldObject.mat.v4[0][3] = 0;
	worldObject.mat.v4[1][3] = 0;
	worldObject.mat.v4[2][3] = 0;
	worldObject.mat.v4[3][0] = 0;
	worldObject.mat.v4[3][1] = 0;
	worldObject.mat.v4[3][2] = 0;
	worldObject.mat.v4[3][3] = 1;

	*went_out_of_bounds = true;
	for(depth = 1; chemin[depth-1].weight > 0 && depth <= MAX_DEPTH; depth ++) {
		rrFactor = 1.0;
		if(depth > MAX_DEPTH_BEFORE_ROULETTE) {
			rrStopProbability = .3;
            if(randFloat(ctx) > rrStopProbability) {
				*went_out_of_bounds = false;
				break;
			}
			rrFactor = rrStopProbability;
		}

		if(!traverseScene(ray, ctx))
			break;

		interaction = ctx->interaction;
		if(vec3_dot(ray.direction, ctx->interaction.n) < 0) {
			interaction.n *= -1;
			interaction.dpdu *= -1;
		}

		// Prepare the matrix
		worldObject.mat.v4[0][0] = interaction.dpdu.x;
		worldObject.mat.v4[0][1] = interaction.dpdv.x;
		worldObject.mat.v4[0][2] = interaction.n.x;

		worldObject.mat.v4[1][0] = interaction.dpdu.y;
		worldObject.mat.v4[1][1] = interaction.dpdv.y;
		worldObject.mat.v4[1][2] = interaction.n.y;

		worldObject.mat.v4[2][0] = interaction.dpdu.z;
		worldObject.mat.v4[2][1] = interaction.dpdv.z;
		worldObject.mat.v4[2][2] = interaction.n.z;

		/*
			On profite du fait que la base (dpdu, dpdv, n) soit une
			base orthonormée pour optimiser le code
		*/
		//matrix_inverse(&worldObject.matInv, worldObject.mat);
		matrix_transpose(&worldObject.matInv, worldObject.mat);

		// Compute the ray's next direction
		wi = vec3_normalize(transform_apply_vector(ray.direction, transform_inverse(worldObject)));
		wo = material_sample_wo(&interaction.m, wi, ctx, current_IOR, &pdf);
		
		if(wo.x == 0 && wo.y == 0 && wo.z == 0)
			break;

		chemin[depth] = (PathElement){
			.pdf = pdf * rrFactor,
			.p = interaction.p - 2 * RAY_EPSILON * interaction.n,
			.weight = material_compute_brdf(&interaction.m, wi, wo, current_IOR, &new_IOR),
			.wi = wi,
			.wo = wo,
			.normalReversed = vec3_dot(ray.direction, ctx->interaction.n) < 0,
            .triangle_id = interaction.triangle_id,
			.worldObject = worldObject,
			.current_IOR = current_IOR
		};

		// Calcul du rayon de la prochaine étape
		ray.origin = chemin[depth].p;
		ray.max_t = MAX_T;
		ray.min_t = RAY_EPSILON;
		ray.direction = vec3_normalize(transform_apply_vector(wo, worldObject));
		current_IOR = new_IOR;
	}

	return depth - 1;
}

Ray genere_rayon_lumineux(Context *ctx, int light_id, Float *pdf)
{
	Ray output;

	output.max_t = MAX_T;
	output.min_t = RAY_EPSILON;

	output.origin = ctx->lights[light_id].pos;

	switch(ctx->lights[light_id].type) {
	case LIGHT_POINT:
		*pdf = .25 * M_1_PI;
		output.direction = UniformSphereSampling(ctx);
		break;

	case LIGHT_AREA:
		output.origin +=
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_width ) *
				ctx->lights[light_id].area_tan +
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_height) *
				ctx->lights[light_id].area_bitan;

		vec3 wo = CosineSampleHemisphere(ctx);
		
		*pdf = M_1_PI * wo.z / (ctx->lights[light_id].area_width *
			ctx->lights[light_id].area_height);

		output.direction =
			wo.x * ctx->lights[light_id].area_tan +
			wo.y * ctx->lights[light_id].area_bitan +
			wo.z * ctx->lights[light_id].area_n;
		break;
	}

	return output;
}

kernel void compute_sample(
	__constant Material *materials,
	__constant vec3 *vertices,
	__constant vec2 *uv,
	__constant Triangle *tris,
	__constant BVHNode *bvh_nodes,
	__constant HostContext *host_ctx,
	__constant Light *lights,
	__global TaskOutput *echantillons
)
{
	int i, j;
	const int id = get_global_id(0);

	Context ctx = {
		.interaction = (SurfaceInteraction) {0},
		.materials = materials,
		.vertices = vertices,
		.uv = uv,
		.tris = tris,
		.bvh_nodes = bvh_nodes,
		.host_ctx = host_ctx,
		.lights = lights
	};

	tinymt32_init(&ctx.twister, host_ctx->seed + id);

	Float pdf;
	PathElement chemin[MAX_DEPTH + 1];

	// L'ID de notre point de départ
	int light_id = randFloat(&ctx) * host_ctx->scene_meta.lights_count;	
	bool went_out_of_bounds;
	int chemin_length = calcul_chemin(&ctx, genere_rayon_lumineux(&ctx, light_id, &pdf), chemin, &went_out_of_bounds);

	pdf *= 1. / (Float) host_ctx->scene_meta.lights_count;
	pdf *= chemin[0].pdf;

	Float weight = chemin[0].weight;
	Float I = lights[light_id].I;

	echantillons[id].went_out_of_bounds = went_out_of_bounds;
	echantillons[id].length = chemin_length;
	for(i = 0; i < chemin_length; i ++) {
		if(i == 0)
			I /= vec3_norm_2(chemin[i+1].p - chemin[i].p);

		weight *= chemin[i+1].weight * fabs(chemin[i+1].wo.z);
		pdf *= chemin[i+1].pdf;

		echantillons[id].triangles_associated[i] = chemin[i+1].triangle_id;
		echantillons[id].values[i] = I * weight / pdf;
	}
}
