#define MAX_DEPTH_BEFORE_ROULETTE	0
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
	__constant unsigned char *maps;
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
	Float weight;
    int triangle_id;
	bool normalReversed;
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

Ray getRayFromLightSource(Context *ctx, int light_id) {
	vec3 w;
	Ray output = (Ray) {
		.origin = ctx->lights[light_id].pos,
		.min_t = RAY_EPSILON,
		.max_t = 1e10
	};

	switch(ctx->lights[light_id].type) {
	case LIGHT_POINT:
		output.direction = UniformSphereSampling(ctx);
		break;

	case LIGHT_AREA:
		w = UniformSampleHemisphere(ctx);
		output.direction = 
			w.x * ctx->lights[light_id].area_tan +
			w.y * ctx->lights[light_id].area_bitan;
		
		if(randFloat(ctx) <= 0.5)
			output.direction += w.z * ctx->lights[light_id].area_n;
		else
			output.direction -= w.z * ctx->lights[light_id].area_n;
		
		output.origin +=
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_width ) * ctx->lights[light_id].area_tan +
			((randFloat(ctx) - .5f) * ctx->lights[light_id].area_height) * ctx->lights[light_id].area_bitan;
		break;
	}

	return output;
}


Float material_computeBRDF(Material *m, vec3 wi, vec3 wo, int lambda)
{
    vec3 valid_wo;
	Float output;
    Float A, B, sigma2, max_cos;
    Float sin_alpha, tan_beta;
    Float sin_theta_wi, sin_theta_wo;
    Float cos_theta_wi, cos_theta_wo;
    Float sin_phi_wi, sin_phi_wo;
    Float cos_phi_wi, cos_phi_wo;

    output = 0;

	cos_theta_wi = wi.z;
	cos_theta_wo = wo.z;
	sin_theta_wi = sqrt(1 - cos_theta_wi * cos_theta_wi);
	sin_theta_wo = sqrt(1 - cos_theta_wo * cos_theta_wo);

    switch(m->type) {
        case MATERIAL_NONE: break;
        case MATERIAL_LAMBERTIAN:
            output = m->l.R * M_1_PI;
            break;
        case MATERIAL_MIRROR:
            output = 1;
            break;
		
        case MATERIAL_GLASS:
			output = 1;
			break;
    }

    return output;
}

vec3 material_sampleWo(Material *m, vec3 wi, Context *ctx)
{
	vec3 wo;
	Float R0, R_theta;
    Float cos_theta_wi;
	Float theta_incident;
    Float theta_refracted;
    Float sin_theta_refracted;

	wo = VEC3(0, 0, 0);
    switch(m->type) {
        case MATERIAL_NONE:
            break;
        case MATERIAL_GLASS:
			wi.z = fabs(wi.z);
            cos_theta_wi = wi.z;
			theta_incident = acos(cos_theta_wi);
            sin_theta_refracted = 1. / m->g.IOR * sin(theta_incident);
			wo = wi;

			R0 = (1. - m->g.IOR) / (1. + m->g.IOR);
			R0 = R0 * R0;
			R_theta = R0 + (1-R0) * (1. - wi.z) * (1. - wi.z) * (1. - wi.z) * (1. - wi.z) * (1. - wi.z);

#if PRINT_DEBUG>=3
			if(ctx->selectedByCursor) {
				printf("Wi.z: " FLOAT_FMT "; Theta_i: " FLOAT_FMT "; sin(Theta_o): " FLOAT_FMT "; IOR: " FLOAT_FMT "; R(theta): " FLOAT_FMT "\n",
					wi.z, theta_incident, sin_theta_refracted, m->g.IOR, R_theta);
			}
#endif	
			if(randFloat(ctx) > R_theta && sin_theta_refracted <= 1 && sin_theta_refracted >= -1) {
				theta_refracted = asin(sin_theta_refracted);
				wo.z -= cos_theta_wi;
				if(cos_theta_wi > 0)
					wo.z += cos(theta_refracted);
				else
					wo.z -= cos(theta_refracted);
				wo = vec3_normalize(wo);
			} else
            	wo.z *= -1;
			break;

        case MATERIAL_MIRROR:
			wo = wi;
            wo.z *= -1;
			break;
        case MATERIAL_LAMBERTIAN:
			wo = CosineSampleHemisphere(ctx);
			wo.z *= -1;
			break;
    }
	return wo;
}

Float light_pdf(Context *ctx, int light_id) {
	switch(ctx->lights[light_id].type) {
	case LIGHT_POINT:
		return .25 * M_1_PI;
	case LIGHT_AREA:
		return .25 * M_1_PI / (ctx->lights[light_id].area_width * ctx->lights[light_id].area_height);
	}
}

Float material_pdf_Wo(Material *m, vec3 wi, vec3 wo, Context *ctx)
{
	switch(m->type) {
        case MATERIAL_NONE:
        case MATERIAL_MIRROR:
        case MATERIAL_GLASS:
            return 1;
        case MATERIAL_LAMBERTIAN:
			return fabs(wo.z);
    }
}

bool interactionRayBBox(BBox3 b, Ray r)
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

	return tmax >= tmin && (tmin > 0 ? tmin : tmax) >= 0 && tmin <= r.max_t;
}

bool interactionRaySphere(Ray r, vec3 center, Float radius)
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

	return t0 > r.min_t && t0 < r.max_t;
}

bool interactionRayLight(Ray r, Light l)
{
	BBox3 b;
	vec3 tmp;
	switch(l.type) {
	case LIGHT_POINT:
		return interactionRaySphere(r, l.pos, 1e-1);
	
	case LIGHT_AREA:
		tmp =
			(Float) (.5)   * l.area_width   * l.area_tan   +
			(Float) (.5)   * l.area_height  * l.area_bitan +
			(Float) (RAY_EPSILON) * l.area_n;
		b = (BBox3) {
			.min = l.pos - tmp,
			.max = l.pos + tmp
		};
		return interactionRayBBox(b, r);
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
	if(t <= r.min_t || t >= r.max_t)
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
		else if(!interactionRayBBox(ctx->bvh_nodes[stack[last_node]].bounds, r))
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

	if(ctx->interaction.m.normal_map.width > 0 && ctx->interaction.m.normal_map.height > 0) {
		i = ctx->interaction.uv.x * ctx->interaction.m.normal_map.width;
		j = ctx->interaction.uv.y * ctx->interaction.m.normal_map.height;

		x = (Float) ctx->maps[ctx->interaction.m.normal_map.start + 3 * (ctx->interaction.m.normal_map.width * j + i) + 0] / 255.;
		y = (Float) ctx->maps[ctx->interaction.m.normal_map.start + 3 * (ctx->interaction.m.normal_map.width * j + i) + 1] / 255.;
		z = (Float) ctx->maps[ctx->interaction.m.normal_map.start + 3 * (ctx->interaction.m.normal_map.width * j + i) + 2] / 255.;
		
		n = x * ctx->interaction.dpdu + y * ctx->interaction.dpdv + z * ctx->interaction.n;
		n = vec3_normalize(n);
	} else
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
		printf("Normal map: %d %d\n",
			ctx->interaction.m.normal_map.width, ctx->interaction.m.normal_map.height);
		printf("Time : " FLOAT_FMT "\n",
			ctx->interaction.time);
	}
#endif

	if(ctx->selected_triangle == ctx->interaction.triangle_id)
		return 1.;
	else
		return fabs(vec3_dot(ray.direction, n));
}


int compute_ray(Ray ray, Context *ctx, PathElement *path, Float w0, int max_depth_before_roulette)
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

    weight = w0;

	path[0] = (PathElement){
		.p = ray.origin,
		.weight = w0,
		.normalReversed = false,
        .triangle_id = -1
	};

	addLastPoint = true;
	for(depth = 0; weight > 0 && depth < MAX_DEPTH; depth ++) {
		rrFactor = 1.0;
		if(depth > max_depth_before_roulette) {
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
			addLastPoint = false;
			break;
		}

		interaction = ctx->interaction;
		if(vec3_dot(ray.direction, ctx->interaction.n) < 0) {
			interaction.n *= -1;
			interaction.dpdu *= -1;
		}

		interaction.n    = vec3_normalize(interaction.n);
		interaction.dpdu = vec3_normalize(interaction.dpdu);
		interaction.dpdv = vec3_normalize(interaction.dpdv);
/*
		if(interaction.m.normal_map.width > 0 && interaction.m.normal_map.height > 0) {
			i = interaction.uv.x * interaction.m.normal_map.width;
			j = interaction.uv.y * interaction.m.normal_map.height;

			x = (Float) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 0] / 255.;
			y = (Float) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 1] / 255.;
			z = (Float) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 2] / 255.;
			interaction.n    = x * ctx->interaction.dpdu + y * ctx->interaction.dpdv + z * ctx->interaction.n;
			interaction.dpdu += interaction.n - ctx->interaction.n;
			interaction.dpdv += interaction.n - ctx->interaction.n;

			interaction.n = vec3_normalize(interaction.n);
			interaction.dpdu = vec3_normalize(interaction.n);
			interaction.dpdv = vec3_normalize(interaction.n);
		}
*/

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
		wo = material_sampleWo(&interaction.m, wi, ctx);
		
		if(wo.x == 0 && wo.y == 0 && wo.z == 0) {
#if PRINT_DEBUG>=3
			if(ctx->selectedByCursor)
				printf("No more WO !\n");
#endif
			break;
		}
		
		weight = weight
			* fabs(wo.z)
			* material_computeBRDF(&interaction.m, wi, wo, ctx->host_ctx->lambda);
			/ material_pdf_Wo(&interaction.m, wi, wo, ctx)
			;

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
			printf("BRDF: " FLOAT_FMT "\n", material_computeBRDF(&interaction.m, wi, wo, ctx->host_ctx->lambda));
			printf("Time : " FLOAT_FMT "\nWeight : " FLOAT_FMT "\nTriangle : %d\n",
				interaction.time, weight, interaction.triangle_id);
		}
#endif

		path[depth+1] = (PathElement){
			.p = interaction.p,
			.weight = weight,
			.normalReversed = vec3_dot(ray.direction, ctx->interaction.n) < 0,
            .triangle_id = interaction.triangle_id
		};

		// Calcul du rayon de la prochaine étape
		ray.origin = interaction.p;
		ray.max_t = ctx->host_ctx->camera.max_t;
		ray.direction = vec3_normalize(transform_apply_vector(wo, worldObject));
	}

	return depth+1;
}

kernel void compute_pixel(
	__constant Material *materials,
	__constant vec3 *vertices,
	__constant vec2 *uv,
	__constant Triangle *tris,
	__constant BVHNode *bvh_nodes,
	__constant HostContext *host_ctx,
	__global int *selected_triangle,
	__constant unsigned char *maps,
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
		.maps = maps,
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
			if(interactionRayLight(ray, ctx.lights[i]))
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
		int cam_path_length = 0;
		int light_path_length = 0;
		PathElement cam_path[MAX_DEPTH];
		PathElement light_path[MAX_DEPTH];
		
#if PRINT_DEBUG>=2
		if(ctx.selectedByCursor)
			printf("Computing camera's path\n");
#endif
		cam_path_length = compute_ray(ray, &ctx, cam_path, 1, 1);


#if PRINT_DEBUG>=2
		if(ctx.selectedByCursor)
			printf("Computing light's path\n");
#endif
		int light_id = randFloat(&ctx) * host_ctx->scene_meta.lights_count;
		light_path_length = compute_ray(getRayFromLightSource(&ctx, light_id), &ctx, light_path, light_pdf(&ctx, light_id) / host_ctx->scene_meta.lights_count, 1);

		Float tmp, totalVal;

		int depth;
		int samplingCount = 0;
		totalVal = 0;
        vec3 dist, wo, wi;
		for(i = 1; i < cam_path_length; i ++) {
			for(j = 0; j < light_path_length; j ++) {
				depth = i + j - 2;
				if((i == 1 && j == 1) || depth < 0 || depth > MAX_DEPTH)
					continue;

			    samplingCount ++;

                dist = light_path[j].p - cam_path[i].p;
                wo = vec3_normalize(dist);
				ray = (Ray) {
					.origin = cam_path[i].p,
					.direction = wo,
					.max_t = vec3_norm(dist),       //On ne va pas plus loin que les 2 points
                    .min_t = RAY_EPSILON            		//NOTE: Purement arbitraire
				};

                wi = vec3_normalize(cam_path[i].p - cam_path[i-1].p);
                
                // TODO: wi et wo dans le plan normal

				tmp = -1;
                Material mat = ctx.materials[ctx.tris[cam_path[i].triangle_id].material];
				if(!traverseScene(ray, &ctx))
					tmp = light_path[j].weight *
						cam_path[i].weight *
                        material_computeBRDF(&mat, wi, wo, ctx.host_ctx->lambda) *                        
						lights[light_id].I;
#if PRINT_DEBUG>=3
				if(ctx.selectedByCursor) {
					printf("Output between %d (" FLOAT_FMT ";" FLOAT_FMT ";" FLOAT_FMT ") and %d (" FLOAT_FMT ";" FLOAT_FMT ";" FLOAT_FMT "): " FLOAT_FMT ";",
                           i, cam_path[i].p.x, cam_path[i].p.y, cam_path[i].p.z,
                           j, light_path[j].p.x, light_path[j].p.y, light_path[j].p.z,
                           tmp);
                    if(tmp == -1)
                        printf("Obstructed by triangle %d, time: " FLOAT_FMT "\n", ctx.interaction.triangle_id, ctx.interaction.time);
                    else
                        printf("\n");
                }
#endif
				if(tmp > 0) {
					totalVal += tmp;
                }
			}
		}
		if(host_ctx->sampleCount == 1)
			raw[index] = 0;
		if(samplingCount > 0)
			raw[index] += totalVal / samplingCount;

#if PRINT_DEBUG>=1
		if(ctx.selectedByCursor)
			printf("Mean weight : " FLOAT_FMT "; Light's path's length: %d; Camera's path's length: %d; Selected triangle: %d\n",
				raw[index] / host_ctx->sampleCount, cam_path_length, light_path_length, ctx.selected_triangle);
#endif

		tmp = raw[index] / host_ctx->sampleCount;
		if(tmp > host_ctx->max_threshold)
			write_imagei(output, on_screen_pos, INT4(255, 0, 0, 0));
		else if(tmp < 0)
			write_imagei(output, on_screen_pos, INT4(0, 255, 0, 0));
		else {
			grayscale = 127 * tmp;
			write_imagei(output, on_screen_pos, INT4(grayscale, grayscale, grayscale, 0));
		} 
	}
}

