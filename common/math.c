#ifdef IN_OPENCL
#include <opencl-c-base.h>
#else
#include "main/typedef.h"
#endif

float lerp(float t, float x1, float x2)
{
	return (1 - t) * x1 + t * x2;
}

typedef struct Mat4x4 Mat4x4;
struct Mat4x4
{
    float c[4][4];
};

Mat4x4 m44_create(float f)
{
    Mat4x4 o;

    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            o.c[i][j] = f;

    return o;
}

Mat4x4 m44_identity()
{
    Mat4x4 o = m44_create(0);

    for(int i = 0; i < 4; i ++)
        o.c[i][i] = 1;
    
    return o;
}

Mat4x4 m44_add(Mat4x4 m1, Mat4x4 m2)
{
    Mat4x4 o;

    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            o.c[i][j] = m1.c[i][j] + m2.c[i][j];

    return o;
}

Mat4x4 m44_mult(Mat4x4 m1, Mat4x4 m2)
{
    Mat4x4 o;

    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            for(int k = 0; k < 4; k ++)
                o.c[i][j] = m1.c[i][k] + m2.c[k][j];

    return o;
}

Mat4x4 m44_smult(Mat4x4 m, double a)
{
    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            m.c[i][j] *= a;

    return m;
}

Mat4x4 m44_transpose(Mat4x4 m)
{
    Mat4x4 o;

    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            o.c[i][j] = m.c[j][i];

    return o;
}

/*
Mat4x4 m44_inverse(Mat4x4 m)
{
    Mat4x4 o;

    for(int i = 0; i < 4; i ++)
        for(int j = 0; j < 4; j ++)
            o.c[i][j] = m.c[j][i];

    return o;
}
*/