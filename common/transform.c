#include "transform.h"
#include "vec3.h"
#include "ray.h"
#include "bbox.h"

Transform transform_inverse(Transform t)
{
    return (Transform) {
        .mat = t.matInv,
        .matInv = t.mat
    };
}

vec3 transform_apply_vector(vec3 p, Transform t)
{
    vec3 output;

    output.x =
        t.mat.s[4*0 + 0] * p.x +
        t.mat.s[4*0 + 1] * p.y +
        t.mat.s[4*0 + 2] * p.z +
        t.mat.s[4*0 + 3];
    output.y =
        t.mat.s[4 * 1 + 0] * p.x +
        t.mat.s[4 * 1 + 1] * p.y +
        t.mat.s[4 * 1 + 2] * p.z +
        t.mat.s[4 * 1 + 3];
    output.z =
        t.mat.s[4 * 2 + 0] * p.x +
        t.mat.s[4 * 2 + 1] * p.y +
        t.mat.s[4 * 2 + 2] * p.z +
        t.mat.s[4 * 2 + 3];
    
    return output;
}
