#include "typedef.h"

#include "vec3.h"

vec3 vec3_add(vec3 v1, vec3 v2)
{
	v1.x += v2.x;
	v1.y += v2.y;
	v1.z += v2.z;
	
	return v1;
}

vec3 vec3_diff(vec3 v1, vec3 v2)
{
	v1.x -= v2.x;
	v1.y -= v2.y;
	v1.z -= v2.z;

	return v1;
}

vec3 vec3_mul(vec3 v1, vec3 v2)
{
	v1.x *= v2.x;
	v1.y *= v2.y;
	v1.z *= v2.z;

	return v1;
}

vec3 vec3_div(vec3 v1, vec3 v2)
{
	v1.x /= v2.x;
	v1.y /= v2.y;
	v1.z /= v2.z;

	return v1;
}

double vec3_dot(vec3 v1, vec3 v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec3 vec3_cross(vec3 v1, vec3 v2)
{
	vec3 o;

	o.x = (v1.y * v2.z) - (v1.z * v2.y);
	o.y = (v1.z * v2.x) - (v1.x * v2.z);
	o.z = (v1.x * v2.y) - (v1.y * v2.x);

	return o;
}

vec3 vec3_smul(double a, vec3 v)
{
	v.x *= a;
	v.y *= a;
	v.z *= a;

	return v;
}

vec3 vec3_sadd(double a, vec3 v)
{
	v.x += a;
	v.y += a;
	v.z += a;

	return v;
}

vec3 vec3_sqrt(vec3 v)
{
	v.x = sqrt(v.x);
	v.y = sqrt(v.y);
	v.z = sqrt(v.z);

	return v;
}

double vec3_norm_2(vec3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

double vec3_norm(vec3 v)
{
	return sqrt(vec3_norm_2(v));
}

vec3 vec3_normalize(vec3 v)
{
	double invNorm;
	double norm2 = vec3_norm_2(v);

	if(norm2 == 0)
		v.x = v.y = v.z = 0;
	else {
		invNorm = 1. / sqrt(norm2);
		v.x *= invNorm;
		v.y *= invNorm;
		v.z *= invNorm;
	}

	return v;
}

void vec3_build_coordonate_system(vec3 v1, vec3 *v2, vec3 *v3)
{
	// TODO: Etudier
	if(fabs(v1.x) > fabs(v1.y)) {
		double invLen = 1. / sqrt(v1.x * v1.x + v1.z * v1.z);
		v2->x = -v1.z * invLen;
		v2->y = 0;
		v2->z = v1.x * invLen;
	}
	else {
		double invLen = 1. / sqrt(v1.y * v1.y + v1.z * v1.z);
		v2->x = 0;
		v2->y = v1.z * invLen;
		v2->z = -v1.y * invLen;
	}

	*v3 = vec3_cross(*v2, v1);
}

vec3 vec3_lerp(vec3 v1, vec3 v2, vec3 t)
{
	vec3 one;
	one.x = one.y = one.z = 1;

	return vec3_add(vec3_mul(vec3_diff(one, t), v1), vec3_mul(t, v2));
}

vec3 vec3_min(vec3 v1, vec3 v2)
{
	vec3 o;

	o.x = fmin(v1.x, v2.x);
	o.y = fmin(v1.y, v2.y);
	o.z = fmin(v1.z, v2.z);

	return o;
}

vec3 vec3_max(vec3 v1, vec3 v2)
{
	vec3 o;

	o.x = fmax(v1.x, v2.x);
	o.y = fmax(v1.y, v2.y);
	o.z = fmax(v1.z, v2.z);

	return o;
}

double vec3_dist_2(vec3 p1, vec3 p2)
{
    return vec3_norm_2(vec3_diff(p1, p2));
}

double vec3_dist(vec3 p1, vec3 p2)
{
    return sqrt(vec3_dist_2(p1, p2));
}

int vec3_max_dimension(vec3 p)
{
	if(p.x > p.y) {
		if(p.x > p.z) return 0;
		else return 2;
	} else {
		if(p.y > p.z) return 1;
		else return 2;
	}
}

vec3 vec3_permute(vec3 p, int kx, int ky, int kz)
{
	return VEC3(
		kx == 0 ? p.x : (kx == 1 ? p.y : p.z),
		ky == 0 ? p.x : (ky == 1 ? p.y : p.z),
		kz == 0 ? p.x : (kz == 1 ? p.y : p.z)
	);
}
