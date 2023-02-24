#include "matrix.h"

#ifndef IN_OPENCL

#include <assert.h>
#include <stdio.h>

#else

#define assert(...)
#define putchar(...)

#endif


void matrix_transpose(Matrix4x4 *out, Matrix4x4 m)
{
    for(int col = 0; col < 4; col++)
        for(int row = 0; row < 4; row++)
			out->s[col + 4*row] = m.s[row + 4*col];
}

void matrix_mult(Matrix4x4 *out, Matrix4x4 m1, Matrix4x4 m2)
{
    // TODO: Utiliser un algo en compléxité < O(n^3)
    for(int col = 0; col < 4; col++)
        for(int row = 0; row < 4; row++) {
			out->s[col + 4*row] = 0;
            for(int i = 0; i < 4; i++)
                out->s[col + 4*row] += m1.s[i + 4*row] * m2.s[col + 4*i];
		}
}

void matrix_print(Matrix4x4 m)
{
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++)
            printf(FLOAT_FMT ",", m.s[j + 4*i]);
        printf("\n");
    }
}

void matrix_swap_rows(Matrix4x4 *m, int i, int j)
{
    for(int col = 0; col < 4; col++) {
        Float temp = m->s[col + 4*i];
        m->s[col + 4*i] = m->s[col + 4*j];
        m->s[col + 4*j] = temp;
    }
}

void matrix_mul_row(Matrix4x4 *m, int row, Float coeff)
{
    for(int col = 0; col < 4; col++)
        m->s[col + 4*row] *= coeff;
}

void matrix_add_row(Matrix4x4 *m, int i, int j, Float coeff)
{
    for(int col = 0; col < 4; col++)
        m->s[col + 4*i] += coeff * m->s[col + 4*j];
}

void matrix_inverse(Matrix4x4 *out, Matrix4x4 m)
{
    Float coeff;

    for(int i = 0; i < 4; i++)
    	for(int j = 0; j < 4; j++)
        	out->s[i + 4*j] = (i==j) ? 1:0;

    for(int col = 0; col < 4; col ++) {
        int row = -1;
        Float max = 0.0f;
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
        matrix_mul_row(out, row, 1.0f / coeff);
        matrix_mul_row(&m,  row, 1.0f / coeff);

        matrix_swap_rows(out, col, row);
        matrix_swap_rows(&m,  col, row);
        
        for(int i = 0; i < 4; i++)
            if(i != col) {
                matrix_add_row(out, i, col, -m.s[col + 4*i]);
                matrix_add_row(&m,  i, col, -m.s[col + 4*i]);
            }
        
        //printf("Row: %d; Col: %d\n", row, col);
        //matrix_print(m); putchar('\n'); putchar('\n');
    }
}
