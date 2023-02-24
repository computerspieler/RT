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

Transform transform_transpose(Transform t)
{
	Transform out;
	matrix_transpose(&out.mat, t.mat);
	matrix_transpose(&out.matInv, t.matInv);
    return out;
}

Transform transform_combine(Transform t1, Transform t2)
{
	Transform out;
	matrix_mult(&out.mat, t1.mat, t2.mat);
	matrix_mult(&out.matInv, t1.matInv, t2.matInv);
    return out;
}

Transform transform_translate(vec3 delta)
{
    return (Transform) {
        .mat = {
            .s = {
                1, 0, 0, delta.x,
                0, 1, 0, delta.y,
                0, 0, 1, delta.z,
                0, 0, 0, 1
            }
        },
        .matInv = {
            .s = {
                1, 0, 0, -delta.x,
                0, 1, 0, -delta.y,
                0, 0, 1, -delta.z,
                0, 0, 0, 1
            }
        },
    };
}

Transform transform_scale(vec3 scale)
{
    return (Transform) {
        .mat = {
            .s = {
                scale.x, 0, 0, 0,
                0, scale.y, 0, 0,
                0, 0, scale.z, 0,
                0, 0, 0, 1
            }
        },
        .matInv = {
            .s = {
                1 / scale.x, 0, 0, 0,
                0, 1 / scale.y, 0, 0,
                0, 0, 1 / scale.z, 0,
                0, 0, 0, 1
            }
        },
    };
}

Transform transform_rotate_x(Float angle)
{
    Float c = cos(angle);
    Float s = sin(angle);
	Transform out;

    Matrix4x4 m = {
        .s = {
            1, 0, 0, 0,
            0, c, -s, 0,
            0, s, c, 0,
            0, 0, 0, 1
        }
    };

	out.mat = m;
	matrix_transpose(&out.matInv, m);
    return out;
}

Transform transform_rotate_y(Float angle)
{
    Float c = cos(angle);
    Float s = sin(angle);
	Transform out;

    Matrix4x4 m = {
        .s = {
            c, 0, s, 0,
            0, 1, 0, 0,
            -s, 0, c, 0,
            0, 0, 0, 1
        }
    };

	out.mat = m;
	matrix_transpose(&out.matInv, m);
    return out;
}

Transform transform_rotate_z(Float angle)
{
    Float c = cos(angle);
    Float s = sin(angle);
	Transform out;

    Matrix4x4 m = {
        .s = {
            c, -s, 0, 0,
            s, c, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        }
    };

	out.mat = m;
	matrix_transpose(&out.matInv, m);
    return out;
}

Transform transform_rotate_axis(Float angle, vec3 axis)
{
    vec3 a = vec3_normalize(axis);
    Float c = cos(angle);
    Float s = sin(angle);
	Transform out;

    Matrix4x4 m = {0};

    // First vector
    m.v4[0][0] = a.x * a.x + (1 - a.x * a.x) * c;
    m.v4[1][0] = a.x * a.y * (1 - c) - a.z * s;
    m.v4[2][0] = a.x * a.z * (1 - c) + a.y * s;

    // Second vector
    m.v4[0][1] = a.y * a.x * (1 - c) + a.z * s;
    m.v4[1][1] = a.y * a.y + (1 - a.y * a.y) * c;
    m.v4[2][1] = a.y * a.z * (1 - c) - a.x * s;

    // Third vector
    m.v4[0][2] = a.z * a.x * (1 - c) - a.y * s;
    m.v4[1][2] = a.z * a.y * (1 - c) + a.x * s;
    m.v4[2][2] = a.z * a.z + (1 - a.z * a.z) * c;

	out.mat = m;
	matrix_transpose(&out.matInv, m);
    return out;
}

Transform transform_look_at(vec3 pos, vec3 look, vec3 up)
{
    Matrix4x4 cameraToWorld;
	Transform out;

    cameraToWorld.v4[0][3] = pos.x;
    cameraToWorld.v4[1][3] = pos.y;
    cameraToWorld.v4[2][3] = pos.z;
    cameraToWorld.v4[3][3] = 1;

    vec3 dir = vec3_normalize(vec3_diff(look, pos));
    vec3 right = vec3_cross(vec3_normalize(up), dir);
    vec3 new_up = vec3_cross(dir, right);

    cameraToWorld.v4[0][0] = right.x;
    cameraToWorld.v4[1][0] = right.y;
    cameraToWorld.v4[2][0] = right.z;
    cameraToWorld.v4[3][0] = 0;

    cameraToWorld.v4[0][1] = new_up.x;
    cameraToWorld.v4[1][1] = new_up.y;
    cameraToWorld.v4[2][1] = new_up.z;
    cameraToWorld.v4[3][1] = 0;

    cameraToWorld.v4[0][2] = dir.x;
    cameraToWorld.v4[1][2] = dir.y;
    cameraToWorld.v4[2][2] = dir.z;
    cameraToWorld.v4[3][2] = 0;

	out.mat = cameraToWorld;
	matrix_inverse(&out.matInv, cameraToWorld);
    return out;
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

vec3 transform_apply_point(vec3 p, Transform t)
{
    vec3 output;

    output.x =
        t.mat.v4[0][0] * p.x +
        t.mat.v4[0][1] * p.y +
        t.mat.v4[0][2] * p.z +
        t.mat.v4[0][3];
    output.y =
        t.mat.v4[1][0] * p.x +
        t.mat.v4[1][1] * p.y +
        t.mat.v4[1][2] * p.z +
        t.mat.v4[1][3];
    output.z =
        t.mat.v4[2][0] * p.x +
        t.mat.v4[2][1] * p.y +
        t.mat.v4[2][2] * p.z +
        t.mat.v4[2][3];
    
    Float w =
        t.mat.v4[3][0] * p.x +
        t.mat.v4[3][1] * p.y +
        t.mat.v4[3][2] * p.z +
        t.mat.v4[3][3];

    return vec3_smul(1/w, output);
}

vec3 transform_apply_normal(vec3 p, Transform t)
{
    vec3 output;

    output.x =
        t.matInv.v4[0][0] * p.x +
        t.matInv.v4[0][1] * p.y +
        t.matInv.v4[0][2] * p.z +
        t.matInv.v4[0][3];
    output.y =
        t.matInv.v4[1][0] * p.x +
        t.matInv.v4[1][1] * p.y +
        t.matInv.v4[1][2] * p.z +
        t.matInv.v4[1][3];
    output.z =
        t.matInv.v4[2][0] * p.x +
        t.matInv.v4[2][1] * p.y +
        t.matInv.v4[2][2] * p.z +
        t.matInv.v4[2][3];
    
    return output;
}

Ray transform_apply_ray(Ray r, Transform t)
{
    Ray o = r;
    o.origin = transform_apply_point(r.origin, t);
    o.direction = transform_apply_vector(r.direction, t);

    // TODO: Offset ray origin to edge of error bounds and compute tMax

    return o;
}

BBox3 transform_apply_bbox3(BBox3 b, Transform t)
{
    vec3 nmin = transform_apply_point(b.min, t);
    BBox3 o = {
        .min = nmin,
        .max = nmin
    };

    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.max.x, b.min.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.min.x, b.max.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.min.x, b.min.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.max.x, b.max.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.max.x, b.min.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        VEC3(b.min.x, b.max.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        b.max, t
    ));

    return o;
}

bool transform_swap_handedness(Transform t)
{
    Float det =
        t.mat.v4[0][0] * (t.mat.v4[1][1] * t.mat.v4[2][2] - t.mat.v4[1][2] * t.mat.v4[2][1]) -
        t.mat.v4[0][1] * (t.mat.v4[1][0] * t.mat.v4[2][2] - t.mat.v4[1][2] * t.mat.v4[2][0]) +
        t.mat.v4[0][2] * (t.mat.v4[1][0] * t.mat.v4[2][1] - t.mat.v4[1][1] * t.mat.v4[2][0])
        ;
    
    return det < 0;
}
