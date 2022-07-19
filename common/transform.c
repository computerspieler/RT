#include "transform.h"
#include "double3.h"
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
    return (Transform) {
        .mat = matrix_transpose(t.matInv),
        .matInv = matrix_transpose(t.mat)
    };
}

Transform transform_translate(double3 delta)
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

Transform transform_scale(double3 scale)
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

Transform transform_rotate_x(double angle)
{
    double c = cos(angle);
    double s = cos(angle);

    Matrix4x4 m = {
        .s = {
            1, 0, 0, 0,
            0, c, -s, 0,
            0, s, c, 0,
            0, 0, 0, 1
        }
    };

    return (Transform) {
        .mat = m,
        .matInv = matrix_transpose(m)
    };
}

Transform transform_rotate_y(double angle)
{
    double c = cos(angle);
    double s = cos(angle);

    Matrix4x4 m = {
        .s = {
            c, 0, s, 0,
            0, 1, 0, 0,
            -s, 0, c, 0,
            0, 0, 0, 1
        }
    };

    return (Transform) {
        .mat = m,
        .matInv = matrix_transpose(m)
    };
}

Transform transform_rotate_z(double angle)
{
    double c = cos(angle);
    double s = cos(angle);

    Matrix4x4 m = {
        .s = {
            c, -s, 0, 0,
            s, c, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0
        }
    };

    return (Transform) {
        .mat = m,
        .matInv = matrix_transpose(m)
    };
}

Transform transform_rotate_axis(double angle, double3 axis)
{
    double3 a = double3_normalize(axis);
    double c = cos(angle);
    double s = cos(angle);

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

    return (Transform) {
        .mat = m,
        .matInv = matrix_transpose(m)
    };
}

Transform transform_look_at(double3 pos, double3 look, double3 up)
{
    Matrix4x4 cameraToWorld;

    cameraToWorld.v4[0][3] = pos.x;
    cameraToWorld.v4[1][3] = pos.y;
    cameraToWorld.v4[2][3] = pos.z;
    cameraToWorld.v4[3][3] = 1;

    double3 dir = double3_normalize(double3_diff(look, pos));
    double3 right = double3_cross(double3_normalize(up), dir);
    double3 new_up = double3_cross(dir, right);

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

    return (Transform) {
        .mat = cameraToWorld,
        .matInv = matrix_inverse(cameraToWorld)
    };
}

double3 transform_apply_vector(double3 p, Transform t)
{
    double3 output;

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
    
    return output;
}

double3 transform_apply_point(double3 p, Transform t)
{
    double3 output;

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
    
    double w =
        t.mat.v4[2][0] * p.x +
        t.mat.v4[2][1] * p.y +
        t.mat.v4[2][2] * p.z +
        t.mat.v4[2][3];

    return double3_smul(1/w, output);
}

double3 transform_apply_normal(double3 p, Transform t)
{
    double3 output;

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
    double3 nmin = transform_apply_point(b.min, t);
    BBox3 o = {
        .min = nmin,
        .max = nmin
    };

    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.max.x, b.min.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.min.x, b.max.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.min.x, b.min.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.max.x, b.max.y, b.min.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.max.x, b.min.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        DOUBLE3(b.min.x, b.max.y, b.max.z), t
    ));
    o = bbox_p_union(o, transform_apply_point(
        b.max, t
    ));

    return o;
}

Transform transform_combine(Transform t1, Transform t2)
{
    return (Transform) {
        .mat = matrix_mult(t1.mat, t2.mat),
        .matInv = matrix_mult(t2.matInv, t1.matInv)
    };
}

bool transform_swap_handedness(Transform t)
{
    double det =
        t.mat.v4[0][0] * (t.mat.v4[1][1] * t.mat.v4[2][2] - t.mat.v4[1][2] * t.mat.v4[2][1]) -
        t.mat.v4[0][1] * (t.mat.v4[1][0] * t.mat.v4[2][2] - t.mat.v4[1][2] * t.mat.v4[2][0]) +
        t.mat.v4[0][2] * (t.mat.v4[1][0] * t.mat.v4[2][1] - t.mat.v4[1][1] * t.mat.v4[2][0])
        ;
    
    return det < 0;
}
