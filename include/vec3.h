#ifndef _VEC3_H_
#define _VEC3_H_

#include "typedef.h"

vec3 vec3_add(vec3 v1, vec3 v2);
vec3 vec3_diff(vec3 v1, vec3 v2);
vec3 vec3_mul(vec3 v1, vec3 v2);
vec3 vec3_div(vec3 v1, vec3 v2);
double vec3_dot(vec3 v1, vec3 v2);
double vec3_adot(vec3 v1, vec3 v2);
vec3 vec3_cross(vec3 v1, vec3 v2);
vec3 vec3_smul(double a, vec3 v);
vec3 vec3_sadd(double a, vec3 v);
vec3 vec3_sqrt(vec3 v);
double vec3_norm_2(vec3 v);
double vec3_norm(vec3 v);
vec3 vec3_normalize(vec3 v);
void vec3_build_coordonate_system(vec3 v1, vec3 *v2, vec3 *v3);
vec3 vec3_lerp(vec3 v1, vec3 v2, vec3 t);
vec3 vec3_min(vec3 v1, vec3 v2);
vec3 vec3_max(vec3 v1, vec3 v2);
double vec3_dist_2(vec3 p1, vec3 p2);
double vec3_dist(vec3 p1, vec3 p2);
int vec3_max_dimension(vec3);
vec3 vec3_permute(vec3 p, int kx, int ky, int kz);

#endif
