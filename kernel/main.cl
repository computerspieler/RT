#define MAX_DEPTH	1

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
	__constant double2 *uv;
	__constant Triangle *tris;
	__constant BVHNode *bvh_nodes;
	__constant unsigned char *maps;

	__global HostContext *host_ctx;
	unsigned int state;
	int2 on_screen_pos;
};

double material_computeBRDF(Material *m, vec3 wi, vec3 wo, int lambda);
vec3 material_sampleWo(Material *m, vec3 wi, Context *ctx);

unsigned int randInt(Context *ctx)
{
	unsigned int x = ctx->state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return ctx->state = x;
}

double randDouble(Context *ctx)
{
	return fabs((double) (randInt(ctx)) / (double) (0xFFFFFFFF));
}


// https://www.rorydriscoll.com/2009/01/07/better-sampling/
vec3 CosineSampleHemisphere(Context *ctx)
{
	const double u1 = randDouble(ctx);
	const double u2 = randDouble(ctx);

	const double r = sqrt(u1);
	const double theta = 2 * M_PI * u2;

	vec3 output;

	output.x = cos(theta) * r;
	output.y = sin(theta) * r;
	output.z = sqrt(max(0., 1. - u1));

    return output;
}

double material_computeBRDF(Material *m, vec3 wi, vec3 wo, int lambda)
{
    vec3 valid_wo;
	double output;
    double A, B, sigma2, max_cos;
    double sin_alpha, tan_beta;
    double sin_theta_wi, sin_theta_wo;
    double cos_theta_wi, cos_theta_wo;
    double sin_phi_wi, sin_phi_wo;
    double cos_phi_wi, cos_phi_wo;

    output = 0;

	cos_theta_wi = wi.z;
	cos_theta_wo = wo.z;
	sin_theta_wi = sqrt(1 - cos_theta_wi * cos_theta_wi);
	sin_theta_wo = sqrt(1 - cos_theta_wo * cos_theta_wo);

    switch(m->type) {
        case MATERIAL_NONE: break;
        case MATERIAL_LAMBERTIAN:
            output = m->l.rho * M_1_PI;
            break;
        case MATERIAL_MIRROR:
            output = 1;
            break;
		
        case MATERIAL_GLASS:
			output = 1;
			break;

        case MATERIAL_OREN_NAYAR:
            sigma2 = m->on.sigma * m->on.sigma;
            A = 1. - (sigma2 / (2. * sigma2 + 0.66));
            B = .45 * sigma2 / (sigma2 + 0.09);

            max_cos = 0;
            if(sin_theta_wi != 0 && sin_theta_wo != 0) {
                cos_phi_wi = (sin_theta_wi == 0) ? 1 : clamp(wi.x / sin_theta_wi, -1., 1.);
                cos_phi_wo = (sin_theta_wo == 0) ? 1 : clamp(wo.x / sin_theta_wo, -1., 1.);
                sin_phi_wi = sqrt(1 - cos_phi_wi * cos_phi_wi);
                sin_phi_wo = sqrt(1 - cos_phi_wo * cos_phi_wo);

                // Note : cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                max_cos = fmax(0., cos_phi_wi * cos_phi_wo + sin_phi_wi * sin_phi_wo);
            }

            if(fabs(cos_theta_wi) > fabs(cos_theta_wo)) {
                sin_alpha = sin_theta_wo;
                tan_beta = sin_theta_wi / fabs(cos_theta_wi);
            } else {
                sin_alpha = sin_theta_wi;
                tan_beta = sin_theta_wo / fabs(cos_theta_wo);
            }

            output = m->on.R * M_1_PI * (A + B * max_cos * sin_alpha * tan_beta);
            break;
    }

    return output;
}

vec3 material_sampleWo(Material *m, vec3 wi, Context *ctx)
{
	vec3 wo;
	double R0, R_theta;
    double cos_theta_wi;
	double theta_incident;
    double theta_refracted;
    double sin_theta_refracted;

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

#ifdef PRINT_DEBUG
			if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y) {
				printf("Wi.z: %f; Theta_i: %f; sin(Theta_o): %f; IOR: %f; R(theta): %f\n",
					wi.z, theta_incident, sin_theta_refracted, m->g.IOR, R_theta);
			}
#endif	
			if(randDouble(ctx) > R_theta && sin_theta_refracted <= 1 && sin_theta_refracted >= -1) {
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
        case MATERIAL_OREN_NAYAR:
			wo = CosineSampleHemisphere(ctx);
			if(wi.z > 0)
				wo.z *= -1;
			break;
    }
	return wo;
}

double material_p_Wo(Material *m, vec3 wi, vec3 wo, Context *ctx)
{
	switch(m->type) {
        case MATERIAL_NONE:
        case MATERIAL_MIRROR:
        case MATERIAL_GLASS:
            return 1;
        case MATERIAL_LAMBERTIAN:
        case MATERIAL_OREN_NAYAR:
			return fabs(wo.z) / M_PI;
    }
}

bool interactionRayBBox(Context *ctx, BBox3 b, Ray r)
{
	vec3 t1, t2;
	double tmin, tmax;

	t1 = (b.min - r.origin) / r.direction;
	t2 = (b.max - r.origin) / r.direction;

	tmin = max(
		max(
			(double) min(t1.x, t2.x),
			(double) min(t1.y, t2.y)
		),  (double) min(t1.z, t2.z)
	);

	tmax = min(
		min(
			(double) max(t1.x, t2.x),
			(double) max(t1.y, t2.y)
		),  (double) max(t1.z, t2.z)
	);

	return tmax >= tmin && (tmin > 0 ? tmin : tmax) >= 0 && tmin < r.max_t;
}

bool interactionRayTriangle(uint t_id, Ray r, Context *ctx)
{
	int kx, ky, kz;
	double e0, e1, e2;
	double b0, b1, b2;
	double det, invDet, t;
	vec3 S;
	vec3 p0, p1, p2;
	vec3 pt0, pt1, pt2;
	double2 uv0, uv1, uv2;
	vec3 ray_dir_transformed;

	p0 = ctx->vertices[ctx->tris[t_id].vertices[0]];
	p1 = ctx->vertices[ctx->tris[t_id].vertices[1]];
	p2 = ctx->vertices[ctx->tris[t_id].vertices[2]];

	kz = vec3_max_dimension(fabs(r.direction));
	kx = (kz + 1) % 3;
	ky = (kx + 1) % 3;
	ray_dir_transformed = vec3_permute(r.direction, kx, ky, kz);

	pt0 = vec3_permute(p0 - r.origin, kx, ky, kz);
	pt1 = vec3_permute(p1 - r.origin, kx, ky, kz);
	pt2 = vec3_permute(p2 - r.origin, kx, ky, kz);
	
	S.z = 1. / ray_dir_transformed.z;
	S.x = -ray_dir_transformed.x * S.z;
	S.y = -ray_dir_transformed.y * S.z;

	pt0.x += S.x * pt0.z;
	pt0.y += S.y * pt0.z;
	pt1.x += S.x * pt1.z;
	pt1.y += S.y * pt1.z;
	pt2.x += S.x * pt2.z;
	pt2.y += S.y * pt2.z;

	e0 = pt1.x * pt2.y - pt1.y * pt2.x;
	e1 = pt2.x * pt0.y - pt2.y * pt0.x;
	e2 = pt0.x * pt1.y - pt0.y * pt1.x;

	if((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
		return false;

	det = e0 + e1 + e2;
	if(det == 0)
		return false;
	
	pt0.z *= S.z;
	pt1.z *= S.z;
	pt2.z *= S.z;

	invDet = 1. / det;
	t = (e0 * pt0.z + e1 * pt1.z + e2 * pt2.z) * invDet;
	if(t < r.min_t || t > r.max_t)
		return false;

	b0 = e0 * invDet;
	b1 = e1 * invDet;
	b2 = e2 * invDet;

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
	int stack[MAX_TREE_DEPTH + 1];
	int seen [MAX_TREE_DEPTH + 1];
	bool interact;
	int last_node = 0;

	interact = false;
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
		else if(!interactionRayBBox(ctx, ctx->bvh_nodes[stack[last_node]].bounds, r))
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


double compute_simple_ray(Ray ray, Context *ctx)
{
	SurfaceInteraction interaction;

	if(!traverseScene(ray, ctx)) {
		if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
			ctx->host_ctx->selected_triangle = -1;

#ifdef PRINT_DEBUG
		if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
			printf("Out !\n");
#endif
		return 0;
	}

	interaction = ctx->interaction;
	if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
		ctx->host_ctx->selected_triangle = interaction.triangle_id;

	interaction.n = vec3_normalize(interaction.n);

	if(ctx->host_ctx->selected_triangle == interaction.triangle_id)
		return 1.;
	else
		return fabs(vec3_dot(ray.direction, interaction.n));
}


double compute_ray(Ray ray, Context *ctx)
{
	vec3 wi, wo;
	int depth;
	Transform worldObject;
	SurfaceInteraction interaction;

	int i, j;
	double out;
    double weight;
	double x, y, z;
	double rrFactor;
	double rrStopProbability;

    weight = 1;
	out = 0;

	for(depth = 0; weight > 0; depth ++) {
		rrFactor = 1.0;
		if(depth > MAX_DEPTH) {
			rrStopProbability = 0.9;
			double d = randDouble(ctx);
            if(d <= rrStopProbability) {
#ifdef PRINT_DEBUG
				if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
					printf("Roulette decided: %f\n", d);
#endif
				break;
			}
			rrFactor = 1.0 / (1.0 - rrStopProbability);
		}

#ifdef PRINT_DEBUG
		if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y) {
			printf("Depth: %d\n", depth);
			printf("Wi: %f %f %f\n",
				ray.direction.x, ray.direction.y, ray.direction.z);
			printf("Origin: %f %f %f\n",
				ray.origin.x, ray.origin.y, ray.origin.z);
		}
#endif

		if(!traverseScene(ray, ctx)) {
#ifdef PRINT_DEBUG
			if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
				printf("Out !\n");
#endif
			break;
		}

		ctx->interaction.n = vec3_normalize(ctx->interaction.n);
		ctx->interaction.dpdu = vec3_normalize(ctx->interaction.dpdu);
		ctx->interaction.dpdv = vec3_normalize(ctx->interaction.dpdv);
		interaction = ctx->interaction;

		if(interaction.m.normal_map.width > 0 && interaction.m.normal_map.height > 0) {
			i = interaction.uv.x * interaction.m.normal_map.width;
			j = interaction.uv.y * interaction.m.normal_map.height;

			x = (double) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 0] / 255.;
			y = (double) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 1] / 255.;
			z = (double) ctx->maps[interaction.m.normal_map.start + 3 * (interaction.m.normal_map.width * j + i) + 2] / 255.;
			interaction.n    = x * ctx->interaction.dpdu + y * ctx->interaction.dpdv + z * ctx->interaction.n;
			interaction.dpdu += interaction.n - ctx->interaction.n;
			interaction.dpdv += interaction.n - ctx->interaction.n;

			interaction.n = vec3_normalize(interaction.n);
			interaction.dpdu = vec3_normalize(interaction.n);
			interaction.dpdv = vec3_normalize(interaction.n);
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
		out += weight * interaction.m.Le;

		// Compute the ray's next direction
		wi = vec3_normalize(transform_apply_vector(ray.direction, transform_inverse(worldObject)));
		wo = material_sampleWo(&interaction.m, wi, ctx);

		
		if(wo.x == 0 && wo.y == 0 && wo.z == 0) {
#ifdef PRINT_DEBUG
			if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y)
				printf("No more WO !\n");
#endif
			return 0;
		}
		
		weight *= fabs(wo.z);
		weight *= material_computeBRDF(&interaction.m, wi, wo, ctx->host_ctx->lambda);
		weight /= material_p_Wo(&interaction.m, wi, wo, ctx);

#ifdef PRINT_DEBUG
		if(ctx->on_screen_pos.x == ctx->host_ctx->mouse_pos.x && ctx->on_screen_pos.y == ctx->host_ctx->mouse_pos.y) {
			printf("Normal: %f %f %f\n",
				interaction.n.x, interaction.n.y, interaction.n.z);
			printf("DPDU: %f %f %f\n",
				interaction.dpdu.x, interaction.dpdu.y, interaction.dpdu.z);
			printf("DPDV: %f %f %f\n",
				interaction.dpdv.x, interaction.dpdv.y, interaction.dpdv.z);
			printf("Hit: %f %f %f\n",
				interaction.p.x, interaction.p.y, interaction.p.z);
			printf("Wi: %f %f %f\n",
				wi.x, wi.y, wi.z);
			printf("Wo: %f %f %f\n",
				wo.x, wo.y, wo.z);
			printf("UV: %f %f\n",
				interaction.uv.x, interaction.uv.y);
			printf("BRDF: %f\n", material_computeBRDF(&interaction.m, wi, wo, ctx->host_ctx->lambda));
			printf("Le: %f\n", interaction.m.Le);
			printf("Time : %f\nOutput : %f\nWeight : %f\n",
				interaction.time, out, weight);
		}
#endif

		// Compute the ray's next step
		ray.origin = interaction.p;
		ray.max_t = ctx->host_ctx->camera.max_t;

		ray.direction = vec3_normalize(transform_apply_vector(wo, worldObject));
	}

    return out;
}

kernel void compute_pixel(
	__constant Material *materials,
	__constant vec3 *vertices,
	__constant double2 *uv,
	__constant Triangle *tris,
	__constant BVHNode *bvh_nodes,
	__global HostContext *host_ctx,
	__constant unsigned char *maps,
	__global double *raw,
	__write_only image2d_t output,
	__global double *surface_irradiance
)
{
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
		.on_screen_pos = on_screen_pos,
		.state =
			host_ctx->seed * on_screen_pos.x * on_screen_pos.y + 1
	};

	double2 normalized_screen_pos = (double2) (
		((double) (2 * on_screen_pos.x) / (double) host_ctx->camera.viewport.x) - 1.0f,
		((double) (2 * on_screen_pos.y) / (double) host_ctx->camera.viewport.y) - 1.0f
	);
	double l = tan(host_ctx->camera.fov) * host_ctx->camera.near;
	vec3 tmp = (vec3)(
		normalized_screen_pos.x * l,
		normalized_screen_pos.y * l,
		host_ctx->camera.near
	);
	Ray ray = (Ray){
		.origin = tmp + host_ctx->camera.pos,
		.direction = vec3_normalize(transform_apply_vector(tmp, host_ctx->camera.rotation_transform)),
		.min_t = 1e-5,
		.max_t = host_ctx->camera.max_t
	};
	
	const int index = on_screen_pos.x + host_ctx->camera.viewport.x * on_screen_pos.y;

	int4 output_final;
	double output_spec;

	if(ctx.host_ctx->simpleView)
		output_spec = compute_simple_ray(ray, &ctx);
	else
		output_spec = compute_ray(ray, &ctx);

	double scale;
	if(host_ctx->simpleView) {
		raw[index] = output_spec;
		scale = 127 * output_spec;
	} else {
		if(host_ctx->sampleCount <= 1)
			raw[index] = 0;	
		raw[index] += output_spec;

#ifdef PRINT_DEBUG
		if(ctx.on_screen_pos.x == ctx.host_ctx->mouse_pos.x && ctx.on_screen_pos.y == ctx.host_ctx->mouse_pos.y) {
			printf("Mean weight : %f\n", raw[index] / host_ctx->sampleCount);
		}
#endif
		scale = 127 * raw[index] / host_ctx->sampleCount;
	}

	output_final = INT4(scale, scale, scale, 0);
	write_imagei(output, on_screen_pos, output_final);
}

