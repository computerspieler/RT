#ifndef _DOUBLE3_H_
#define _DOUBLE3_H_

#include "typedef.h"

double3 double3_add(double3 v1, double3 v2);
double3 double3_diff(double3 v1, double3 v2);
double3 double3_mul(double3 v1, double3 v2);
double3 double3_div(double3 v1, double3 v2);
double double3_dot(double3 v1, double3 v2);
double double3_adot(double3 v1, double3 v2);
double3 double3_cross(double3 v1, double3 v2);
double3 double3_smul(double a, double3 v);
double3 double3_sadd(double a, double3 v);
double3 double3_sqrt(double3 v);
double double3_norm_2(double3 v);
double double3_norm(double3 v);
double3 double3_normalize(double3 v);
void double3_build_coordonate_system(double3 v1, double3 *v2, double3 *v3);
double3 double3_lerp(double3 v1, double3 v2, double3 t);
double3 double3_min(double3 v1, double3 v2);
double3 double3_max(double3 v1, double3 v2);
double double3_dist_2(double3 p1, double3 p2);
double double3_dist(double3 p1, double3 p2);

#endif
