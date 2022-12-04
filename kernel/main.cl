#ifndef IN_OPENCL
#define IN_OPENCL

// This part will never be compiled
// This is just for the completion
#define MAX_TREE_DEPTH 0
#include <opencl-c-base.h>

#endif

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
#include <math.h>
#endif

typedef struct Context Context;
struct Context
{
	SurfaceInteraction interaction;
	__constant Material *materials;
	__constant vec3 *vertices;
	__constant vec3 *normals;
	__constant double2 *uv;
	__constant Triangle *tris;
	__constant BVHNode *bvh_nodes;
	__constant HostContext *host_ctx;
	
	unsigned int state;
	int2 on_screen_pos;
};

double material_computeBRDF(Material *m, vec3 wi, vec3 wo, double lambda);
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
	const double phi = 2 * M_PI * u2;

	vec3 output;

	output.x = cos(phi) * r;
	output.y = sin(phi) * r;
	output.z = sqrt(fmax(0., 1 - u1));

    return output;
}
vec3 UniformSphereSampling(Context *ctx)
{
	const double u1 = 2 * M_PI * randDouble(ctx);
	const double u2 = 2 * M_PI * randDouble(ctx);

	vec3 output;

	const double r = cos(u1);

	output.x = cos(u2) * r;
	output.y = sin(u2) * r;
	output.z = sin(u1);

    return output;
}

double light_Power(__constant Light *l)
{
	switch(l->type) {
	case LIGHT_POINT:
		return 4 * M_PI * l->I;
	}
}
double light_LE(__constant Light *l, int lambda)
{
	switch(l->type) {
	case LIGHT_POINT:
		return l->I;
	}
}

double material_computeBRDF(Material *m, vec3 wi, vec3 wo, double lambda)
{
    vec3 valid_wo;
	double output;
	double pdf;
    double A, B, sigma2, max_cos;
    double sin_alpha, tan_beta;
    double sin_theta_wi, sin_theta_wo;
    double cos_theta_wi, cos_theta_wo;
    double sin_phi_wi, sin_phi_wo;
    double cos_phi_wi, cos_phi_wo;

	pdf = 0;
    output = 0;

	cos_theta_wi = wi.z;
	cos_theta_wo = wo.z;
	sin_theta_wi = sqrt(1 - cos_theta_wi * cos_theta_wi);
	sin_theta_wo = sqrt(1 - cos_theta_wo * cos_theta_wo);

    switch(m->type) {
        case MATERIAL_NONE: break;
        case MATERIAL_LAMBERTIAN:
            output = m->l.rho * M_1_PI;
			pdf = cos_theta_wo * M_1_PI;
            break;
        case MATERIAL_MIRROR:
            valid_wo = wi - VEC3(0, 0, 2 * cos_theta_wi);
            if(wo.x == valid_wo.x && wo.y == valid_wo.y && wo.z == valid_wo.z) {
                output = 1;
				pdf = 1;
			}
            break;
		
        case MATERIAL_GLASS:
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
			pdf = cos_theta_wo * M_1_PI;
            break;
    }

    return output / pdf;
}

vec3 material_sampleWo(Material *m, vec3 wi, Context *ctx)
{
    double cos_theta_wi;
    switch(m->type) {
        case MATERIAL_NONE:
            return VEC3(0, 0, 0);
        case MATERIAL_MIRROR:
            cos_theta_wi = wi.z;
            return wi - VEC3(0, 0, 2 * cos_theta_wi);
        case MATERIAL_GLASS:
            return VEC3(0, 0, 0);
        case MATERIAL_LAMBERTIAN:
        case MATERIAL_OREN_NAYAR:
			return CosineSampleHemisphere(ctx);
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

void compute_ray(Ray ray, Context *ctx, Spectrum *out, vec3 *last_pos)
{
	int i;
	vec3 wi, wo;
	double cos_angle;
	int depth;
	Transform worldObject;
	SurfaceInteraction interaction;

	double rrFactor;
	double rrStopProbability;

	for(i = 0; i < SPECTRUM_SIZE; i ++)
		(*out)[i] = 0;

	for(depth = 0; ; depth ++) {
		rrFactor = 1.0;
		if(depth >= MAX_DEPTH) {
			rrStopProbability = .3;
			if(randDouble(ctx) > rrStopProbability)
				break;
			rrFactor = 1.0 / (1.0 - rrStopProbability);
		}

#ifdef PRINT_DEBUG
		if(ctx->on_screen_pos.x == 160 && ctx->on_screen_pos.y == 160) {
			printf("Depth: %d\n", depth);
			printf("Wi: %lf %lf %lf\n",
				ray.direction.x, ray.direction.y, ray.direction.z);
			printf("Origin: %lf %lf %lf\n",
				ray.origin.x, ray.origin.y, ray.origin.z);
		}
#endif

		if(!traverseScene(ray, ctx)) {
#ifdef PRINT_DEBUG
			if(ctx->on_screen_pos.x == 160 && ctx->on_screen_pos.y == 160)
				printf("Out !\n");
#endif
			break;
		}

		interaction = ctx->interaction;
		cos_angle = vec3_dot(interaction.n, ray.direction);
		if(cos_angle > 0) {
			interaction.n *= -1;
			interaction.dpdu *= -1;
		} else
			cos_angle *= -1;
		
		interaction.dpdu = vec3_normalize(interaction.dpdu);
		interaction.dpdv = vec3_normalize(interaction.dpdv);

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

		worldObject.mat.v4[3][3] = 1;

		worldObject.matInv = matrix_inverse(worldObject.mat);

		double r2 = max(1., ctx->interaction.time * ctx->interaction.time);
		for(i = 0; i < SPECTRUM_SIZE; i ++)
			*out[i] *= cos_angle * rrFactor / r2;
		
		// Compte the ray's next direction
		wi = transform_apply_vector(ray.direction, transform_inverse(worldObject));
		wo = material_sampleWo(&interaction.m, wi, ctx);
		
		if(wo.x == 0 && wo.y == 0 && wo.z == 0)
			break;

		for(i = 0; i < SPECTRUM_SIZE; i ++)
			*out[i] *= material_computeBRDF(&interaction.m, wi, wo, i * SPECTRUM_STEP + SPECTRUM_SIZE);

#ifdef PRINT_DEBUG
		if(ctx->on_screen_pos.x == 160 && ctx->on_screen_pos.y == 160) {
			printf("Normal: %lf %lf %lf\n",
				interaction.n.x, interaction.n.y, interaction.n.z);
			printf("DPDU: %lf %lf %lf\n",
				interaction.dpdu.x, interaction.dpdu.y, interaction.dpdu.z);
			printf("DPDV: %lf %lf %lf\n",
				interaction.dpdv.x, interaction.dpdv.y, interaction.dpdv.z);
			printf("Hit: %lf %lf %lf\n",
				interaction.p.x, interaction.p.y, interaction.p.z);
			printf("Wi: %lf %lf %lf\n",
				wi.x, wi.y, wi.z);
			printf("Wo: %lf %lf %lf\n",
				wo.x, wo.y, wo.z);
			printf("cos(theta): %lf\nTime : %lf\n",
				cos_angle, interaction.time);
		}
#endif

		// Compute the ray's next step
		ray.origin = interaction.p;
		ray.max_t = ctx->host_ctx->camera.max_t;

		ray.direction = vec3_normalize(transform_apply_vector(wo, worldObject));
	}

	if(last_pos)
		*last_pos = ray.origin;
}

kernel void compute_pixel(
	__constant Material *materials,
	__constant vec3 *vertices,
	__constant vec3 *normals,
	__constant double2 *uv,
	__constant Triangle *tris,
	__constant BVHNode *bvh_nodes,
	__constant Light *lights,
	__constant HostContext *host_ctx,
	__global double *raw,
	__write_only image2d_t output
)
{
	const int2 on_screen_pos = {get_global_id(0), get_global_id(1)};

	Context ctx = {
		.interaction = (SurfaceInteraction) {0},
		.materials = materials,
		.vertices = vertices,
		.normals = normals,
		.uv = uv,
		.tris = tris,
		.bvh_nodes = bvh_nodes,
		.host_ctx = host_ctx,
		.on_screen_pos = on_screen_pos,
		.state =
			host_ctx->seed * on_screen_pos.x * on_screen_pos.y
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
	
	float4 pixel_raw;
	int4 output_final;

	int state, i, j;
	double weight_max, weight;
	Spectrum cam_weight, light_weight, output_spec;

	vec3 output_color = VEC3(0, 0, 0);
	vec3 output_cam_last_pos;
	vec3 output_light_last_pos;

	const int index = 4 * (on_screen_pos.x + host_ctx->camera.viewport.x * on_screen_pos.y);

	pixel_raw = FLOAT4(0, 0, 0, 0);
	if(host_ctx->sampleCount != 1)
		pixel_raw = FLOAT4(
			raw[index + 0],
			raw[index + 1],
			raw[index + 2],
			raw[index + 3]
		);
	
	state = ctx.state;
	weight_max = 0;
	compute_ray(ray, &ctx, &cam_weight, &output_cam_last_pos);

#ifdef PRINT_DEBUG
	if(ctx.on_screen_pos.x == 160 && ctx.on_screen_pos.y == 160) {
		printf("Cam max weight: %lf\n", weight_max);
	}
#endif

	for(i = 0; i < host_ctx->scene_meta.lights_count; i ++) {
		if(lights[i].I == 0)
			continue;

		ray.origin = lights[i].pos;
		ray.direction = UniformSphereSampling(&ctx);
		compute_ray(ray, &ctx, &light_weight, &output_cam_last_pos);

		ray.direction = output_light_last_pos - output_cam_last_pos;
		ray.max_t = vec3_norm(ray.direction);
		ray.direction = vec3_normalize(ray.direction);
		ray.origin = output_cam_last_pos;
		if(traverseScene(ray, &ctx))
			break;
			
		spectrum_mult(light_weight, light_weight, cam_weight);
		spectrum_add(light_weight, light_weight, lights[i].light);
		spectrum_mult_c(output_spec, output_spec, light_LE(&lights[i], 0));
		spectrum_add(output_spec, output_spec, light_weight);

#ifdef PRINT_DEBUG
		if(ctx.on_screen_pos.x == 160 && ctx.on_screen_pos.y == 160) {
			printf("Light %d's light: %lf\n", i, lights[i].I);
		}
#endif
	}

	weight_max = 0;
	for(i = 0; i < SPECTRUM_SIZE; i ++) {
		output_color += spectrum_to_rgb(i * SPECTRUM_STEP + SPECTRUM_START) * output_spec[i];
		if(output_spec[i] > weight_max)
			weight_max = output_spec[i];
	}

#ifdef PRINT_DEBUG
	if(ctx.on_screen_pos.x == 160 && ctx.on_screen_pos.y == 160) {
		printf("Light max weight: %lf\n", weight_max);
	}
#endif

	output_color /= host_ctx->scene_meta.lights_count;

	pixel_raw += FLOAT4(
		output_color.x, output_color.y, output_color.z, 0
	);

	raw[index + 0] = pixel_raw.x;
	raw[index + 1] = pixel_raw.y;
	raw[index + 2] = pixel_raw.z;
	raw[index + 3] = pixel_raw.w;

	double scale = 255. / (double) host_ctx->sampleCount;
	output_final = INT4(
		pixel_raw.z * scale,
		pixel_raw.y * scale,
		pixel_raw.x * scale,
		0
	);

#ifdef PRINT_DEBUG
	if(on_screen_pos.x == 160 && on_screen_pos.y == 160)
		output_final = INT4(0, 255, 255, 0);
#endif

	write_imagei(output, on_screen_pos, output_final);
} 	
