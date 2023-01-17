#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "typedef.h"


typedef union Matrix4x4 Matrix4x4;
union Matrix4x4
{
    double16 f;
    double s[16];
    double v4[4][4];
};

void matrix_print (Matrix4x4 m);
void matrix_inverse (Matrix4x4 *output, Matrix4x4 m);
void matrix_transpose (Matrix4x4 *output, Matrix4x4 m);
void matrix_mult (Matrix4x4 *output, Matrix4x4 m1, Matrix4x4 m2);

#endif
