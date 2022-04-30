#ifndef _FLOAT3_H_
#define _FLOAT3_H_

#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

float3 float3_add(float3 v1, float3 v2);
float3 float3_diff(float3 v1, float3 v2);
float3 float3_mul(float3 v1, float3 v2);
float3 float3_div(float3 v1, float3 v2);
double float3_dot(float3 v1, float3 v2);
double float3_adot(float3 v1, float3 v2);
float3 float3_cross(float3 v1, float3 v2);
float3 float3_smul(double a, float3 v);
double float3_norm_2(float3 v);
double float3_norm(float3 v);
float3 float3_normalize(float3 v);
void float3_build_coordonate_system(float3 v1, float3 *v2, float3 *v3);
float3 float3_lerp(float3 v1, float3 v2, float3 t);
float3 float3_min(float3 v1, float3 v2);
float3 float3_max(float3 v1, float3 v2);
double float3_dist_2(float3 p1, float3 p2);
double float3_dist(float3 p1, float3 p2);

#endif