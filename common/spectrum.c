#include "spectrum.h"

#ifndef IN_OPENCL

#include <assert.h>
#include <stdio.h>

#else

#define assert(...)
#define putchar(...)

#endif

double spectrum_sum(Spectrum *s)
{
    double output = 0;
    for(int i = 0; i < SPECTRUM_SIZE; i++)
        output += *s[i];
    return output;
}

double spectrum_avg(Spectrum *s)
{
    return spectrum_sum(s) / SPECTRUM_SIZE;
}