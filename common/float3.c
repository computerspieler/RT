#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

#include "float3.h"

float3 float3_add(float3 v1, float3 v2)
{
	v1.x += v2.x;
	v1.y += v2.y;
	v1.z += v2.z;
	
	return v1;
}

float3 float3_diff(float3 v1, float3 v2)
{
	v1.x -= v2.x;
	v1.y -= v2.y;
	v1.z -= v2.z;

	return v1;
}

float3 float3_mul(float3 v1, float3 v2)
{
	v1.x *= v2.x;
	v1.y *= v2.y;
	v1.z *= v2.z;

	return v1;
}

float3 float3_div(float3 v1, float3 v2)
{
	v1.x /= v2.x;
	v1.y /= v2.y;
	v1.z /= v2.z;

	return v1;
}

double float3_dot(float3 v1, float3 v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double float3_adot(float3 v1, float3 v2)
{
	return fabs(float3_dot(v1, v2));
}

float3 float3_cross(float3 v1, float3 v2)
{
	float3 o;

	o.x = v1.y * v2.z - v1.z * v2.y;
	o.y = v1.z * v2.x - v1.x * v2.z;
	o.z = v1.x * v2.y - v1.y * v2.x;

	return o;
}

float3 float3_smul(double a, float3 v)
{
	v.x *= a;
	v.y *= a;
	v.z *= a;

	return v;
}

double float3_norm_2(float3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

double float3_norm(float3 v)
{
	return sqrt(float3_norm_2(v));
}

float3 float3_normalize(float3 v)
{
	double norm = float3_norm(v);

	if(norm == 0) {
		v.x = 1;
		v.y = v.z = 0;
	} else {
		v.x /= norm;
		v.y /= norm;
		v.z /= norm;
	}

	return v;
}

void float3_build_coordonate_system(float3 v1, float3 *v2, float3 *v3)
{
	// TODO: Etudier
	if(fabs(v1.x) > fabs(v1.y)) {
		float invLen = rsqrt(v1.x * v1.x + v1.z * v1.z);
		v2->x = -v1.z * invLen;
		v2->y = 0;
		v2->z = v1.x * invLen;
	}
	else {
		float invLen = rsqrt(v1.y * v1.y + v1.z * v1.z);
		v2->x = 0;
		v2->y = -v1.z * invLen;
		v2->z = -v1.y * invLen;
	}

	*v3 = float3_cross(v1, *v2);
}

float3 float3_lerp(float3 v1, float3 v2, float3 t)
{
	float3 one;
	one.x = one.y = one.z = 1;

	return float3_add(float3_mul(float3_diff(one, t), v1), float3_mul(t, v2));
}

float3 float3_min(float3 v1, float3 v2)
{
	float3 o;

	o.x = fmin(v1.x, v2.x);
	o.y = fmin(v1.y, v2.y);
	o.z = fmin(v1.z, v2.z);

	return o;
}

float3 float3_max(float3 v1, float3 v2)
{
	float3 o;

	o.x = fmax(v1.x, v2.x);
	o.y = fmax(v1.y, v2.y);
	o.z = fmax(v1.z, v2.z);

	return o;
}

double float3_dist_2(float3 p1, float3 p2)
{
    return float3_norm_2(float3_diff(p1, p2));
}

double float3_dist(float3 p1, float3 p2)
{
    return sqrt(float3_dist_2(p1, p2));
}
