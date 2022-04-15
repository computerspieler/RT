#include <math.h>

#include "vec3.h"

Vec3 vec3_create(double x, double y, double z)
{
	return (Vec3) {
		.x = x,
		.y = y,
		.z = z
	};
}

Vec3 vec3_add(Vec3 v1, Vec3 v2)
{
	return (Vec3) {
		.x = v1.x + v2.x,
		.y = v1.y + v2.y,
		.z = v1.z + v2.z
	};
}

Vec3 vec3_diff(Vec3 v1, Vec3 v2)
{
	return (Vec3) {
		.x = v1.x - v2.x,
		.y = v1.y - v2.y,
		.z = v1.z - v2.z
	};
}

double vec3_dot(Vec3 v1, Vec3 v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double vec3_adot(Vec3 v1, Vec3 v2)
{
	return fabs(vec3_dot(v1, v2));
}

Vec3 vec3_cross(Vec3 v1, Vec3 v2)
{
	return (Vec3) {
		.x = v1.y * v2.z - v1.z * v2.y,
		.y = v1.z * v2.x - v1.x * v2.z,
		.z = v1.x * v2.y - v1.y * v2.x,
	};
}

Vec3 vec3_smul(double a, Vec3 v)
{
	return (Vec3) {
		.x = v.x * a,
		.y = v.y * a,
		.z = v.z * a
	};
}

double vec3_norm_2(Vec3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

double vec3_norm(Vec3 v)
{
	return sqrt(vec3_norm_2(v));
}

Vec3 vec3_normalize(Vec3 v)
{
	double norm = vec3_norm(v);

	if(norm == 0)
		return (Vec3) {1, 0, 0};

	return (Vec3) {
		.x = v.x / norm,
		.y = v.y / norm,
		.z = v.z / norm
	};
}

void vec3_build_coordonate_system(Vec3 v1, Vec3 *v2, Vec3 *v3)
{
	// TODO: Etudier
	if(fabs(v1.x) > fabs(v1.y)) {
		double invLen = 1.0 / sqrtf(v1.x * v1.x + v1.z * v1.z);
		*v2 = (Vec3) {.x = -v1.z * invLen, .y = 0, .z = v1.x * invLen};
	}
	else {
		double invLen = 1.0 / sqrtf(v1.y * v1.y + v1.z * v1.z);
		*v2 = (Vec3) {.x = 0, .y = -v1.z * invLen, .z = -v1.y * invLen};
	}

	*v3 = vec3_cross(v1, *v2);
}

Vec3 vec3_lerp(Vec3 v1, Vec3 v2, double t)
{
	return vec3_add(vec3_smul((1 - t), v1), vec3_smul(t, v2));
}

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

Vec3 vec3_min(Vec3 v1, Vec3 v2)
{
	return (Vec3) {
		.x = min(v1.x, v2.x),
		.y = min(v1.y, v2.y),
		.z = min(v1.z, v2.z)
	};
}

Vec3 vec3_max(Vec3 v1, Vec3 v2)
{
	return (Vec3) {
		.x = max(v1.x, v2.x),
		.y = max(v1.y, v2.y),
		.z = max(v1.z, v2.z)
	};
}