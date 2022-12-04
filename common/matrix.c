#include "matrix.h"

#ifndef IN_OPENCL

#include <assert.h>
#include <stdio.h>

#else

#define assert(...)
#define putchar(...)

#endif


Matrix4x4 matrix_transpose(Matrix4x4 m)
{
    return (Matrix4x4){
        .s = {
            m.v4[0][0], m.v4[1][0], m.v4[2][0], m.v4[3][0],
            m.v4[0][1], m.v4[1][1], m.v4[2][1], m.v4[3][1],
            m.v4[0][2], m.v4[1][2], m.v4[2][2], m.v4[3][2],
            m.v4[0][3], m.v4[1][3], m.v4[2][3], m.v4[3][3]
        }
    };
}

Matrix4x4 matrix_mult(Matrix4x4 m1, Matrix4x4 m2)
{
    // TODO: Utiliser un algo en compléxité < O(n^3)
    Matrix4x4 output = (Matrix4x4) {0};

    for(int col = 0; col < 4; col++)
        for(int row = 0; row < 4; row++)
            for(int i = 0; i < 4; i++)
                output.s[col + 4*row] += m1.s[i + 4*row] * m2.s[col + 4*i];
    
    return output;
}

void matrix_print(Matrix4x4 m)
{
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++)
            printf("%f,", m.s[j + 4*i]);
        printf("\n");
    }
}

void matrix_swap_rows(Matrix4x4 *m, int i, int j)
{
    for(int col = 0; col < 4; col++) {
        double temp = m->s[col + 4*i];
        m->s[col + 4*i] = m->s[col + 4*j];
        m->s[col + 4*j] = temp;
    }
}

void matrix_mul_row(Matrix4x4 *m, int row, double coeff)
{
    for(int col = 0; col < 4; col++)
        m->s[col + 4*row] *= coeff;
}

void matrix_add_row(Matrix4x4 *m, int i, int j, double coeff)
{
    for(int col = 0; col < 4; col++)
        m->s[col + 4*i] += coeff * m->s[col + 4*j];
}

Matrix4x4 matrix_inverse(Matrix4x4 m)
{
    double coeff;
    Matrix4x4 output = (Matrix4x4) {0};

    for(int i = 0; i < 4; i++)
        output.s[i + 4*i] = 1;

    for(int col = 0; col < 4; col ++) {
        int row = -1;
        double max = 0.0f;
        for(int i = col; i < 4; i ++) {
            if(m.s[col + 4 * i] == 0.0f)
                continue;

            if(row == -1 || m.s[col + 4 * i] > max) {
                row = i;
                max = m.s[col + 4 * i];
            }
        }
        
        assert(row != -1);
        
        coeff = m.s[col + 4*row];
        matrix_mul_row(&output, row, 1.0f / coeff);
        matrix_mul_row(&m,      row, 1.0f / coeff);

        matrix_swap_rows(&output, col, row);
        matrix_swap_rows(&m,      col, row);
        
        for(int i = 0; i < 4; i++)
            if(i != col) {
                matrix_add_row(&output, i, col, -m.s[col + 4*i]);
                matrix_add_row(&m,      i, col, -m.s[col + 4*i]);
            }
        

        //printf("Row: %d; Col: %d\n", row, col);
        //matrix_print(m); putchar('\n'); putchar('\n');
    }

    return output;
}
