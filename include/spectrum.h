#ifndef _SPECTRUM_H_
#define _SPECTRUM_H_

#include "typedef.h"

#define SPECTRUM_START  (350)
#define SPECTRUM_END    (700)
#define SPECTRUM_STEP   (80)
#define SPECTRUM_SIZE   ((SPECTRUM_END-SPECTRUM_START) / SPECTRUM_STEP)

typedef double Spectrum[SPECTRUM_SIZE];

#define spectrum_add(sout, s1, s2) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] + s2[i];

#define spectrum_sub(sout, s1, s2) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] - s2[i];

#define spectrum_div(sout, s1, s2) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] / s2[i];

#define spectrum_mult(sout, s1, s2) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] * s2[i];

#define spectrum_add_c(sout, s1, c) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] + c;

#define spectrum_sub_c(sout, s1, c) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] - c;

#define spectrum_div_c(sout, s1, c) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] / c;

#define spectrum_mult_c(sout, s1, c) \
    for(int i = 0; i < SPECTRUM_SIZE; i ++) \
        sout[i] = s1[i] * c;

#define spectrum_val(s, lambda) \
    s[(lambda - SPECTRUM_START) / SPECTRUM_STEP]

double spectrum_sum(Spectrum *s);
double spectrum_avg(Spectrum *s);

double3 spectrum_to_xyz(int lambda);
double3 spectrum_to_rgb(int lambda);

#endif
