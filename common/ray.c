#include "typedef.h"

#include "double3.h"
#include "ray.h"

RayDifferential ray_scale_differentials(RayDifferential r, double s)
{
    return (RayDifferential) {
        .ray = r.ray,
        .rx_origin = double3_add(r.ray.origin, double3_smul(s, double3_diff(r.rx_origin, r.ray.origin))),
        .ry_origin = double3_add(r.ray.origin, double3_smul(s, double3_diff(r.ry_origin, r.ray.origin))),
        .rx_direction = double3_add(r.ray.direction, double3_smul(s, double3_diff(r.rx_direction, r.ray.direction))),
        .ry_direction = double3_add(r.ray.direction, double3_smul(s, double3_diff(r.ry_direction, r.ray.direction)))
    };
}
