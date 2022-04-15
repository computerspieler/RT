#include <math.h>

#include "point.h"

double point_dist_2(Point p1, Point p2)
{
    return vec3_norm_2(vec3_diff(p1, p2));
}

double point_dist(Point p1, Point p2)
{
    return sqrt(point_dist_2(p1, p2));
}