#include "typedef.h"

#include "double3.h"

double3 double3_add(double3 v1, double3 v2)
{
	v1.x += v2.x;
	v1.y += v2.y;
	v1.z += v2.z;
	
	return v1;
}

double3 double3_diff(double3 v1, double3 v2)
{
	v1.x -= v2.x;
	v1.y -= v2.y;
	v1.z -= v2.z;

	return v1;
}

double3 double3_mul(double3 v1, double3 v2)
{
	v1.x *= v2.x;
	v1.y *= v2.y;
	v1.z *= v2.z;

	return v1;
}

double3 double3_div(double3 v1, double3 v2)
{
	v1.x /= v2.x;
	v1.y /= v2.y;
	v1.z /= v2.z;

	return v1;
}

double double3_dot(double3 v1, double3 v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double double3_adot(double3 v1, double3 v2)
{
	return fabs(double3_dot(v1, v2));
}

double3 double3_cross(double3 v1, double3 v2)
{
	double3 o;

	o.x = v1.y * v2.z - v1.z * v2.y;
	o.y = v1.z * v2.x - v1.x * v2.z;
	o.z = v1.x * v2.y - v1.y * v2.x;

	return o;
}

double3 double3_smul(double a, double3 v)
{
	v.x *= a;
	v.y *= a;
	v.z *= a;

	return v;
}

double3 double3_sadd(double a, double3 v)
{
	v.x += a;
	v.y += a;
	v.z += a;

	return v;
}

double3 double3_sqrt(double3 v)
{
	v.x = sqrt(v.x);
	v.y = sqrt(v.y);
	v.z = sqrt(v.z);

	return v;
}

double double3_norm_2(double3 v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

double double3_norm(double3 v)
{
	return sqrt(double3_norm_2(v));
}

double3 double3_normalize(double3 v)
{
	double norm = double3_norm(v);

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

void double3_build_coordonate_system(double3 v1, double3 *v2, double3 *v3)
{
	// TODO: Etudier
	if(fabs(v1.x) > fabs(v1.y)) {
		double invLen = rsqrt(v1.x * v1.x + v1.z * v1.z);
		v2->x = -v1.z * invLen;
		v2->y = 0;
		v2->z = v1.x * invLen;
	}
	else {
		double invLen = rsqrt(v1.y * v1.y + v1.z * v1.z);
		v2->x = 0;
		v2->y = -v1.z * invLen;
		v2->z = -v1.y * invLen;
	}

	*v3 = double3_cross(v1, *v2);
}

double3 double3_lerp(double3 v1, double3 v2, double3 t)
{
	double3 one;
	one.x = one.y = one.z = 1;

	return double3_add(double3_mul(double3_diff(one, t), v1), double3_mul(t, v2));
}

double3 double3_min(double3 v1, double3 v2)
{
	double3 o;

	o.x = fmin(v1.x, v2.x);
	o.y = fmin(v1.y, v2.y);
	o.z = fmin(v1.z, v2.z);

	return o;
}

double3 double3_max(double3 v1, double3 v2)
{
	double3 o;

	o.x = fmax(v1.x, v2.x);
	o.y = fmax(v1.y, v2.y);
	o.z = fmax(v1.z, v2.z);

	return o;
}

double double3_dist_2(double3 p1, double3 p2)
{
    return double3_norm_2(double3_diff(p1, p2));
}

double double3_dist(double3 p1, double3 p2)
{
    return sqrt(double3_dist_2(p1, p2));
}
