#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_

#ifndef IN_OPENCL

#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <CL/opencl.h>

#include <math.h>
#include <stdbool.h>

#define kernel 
#define __global 
#define __write_only 

typedef cl_double16 double16;
typedef cl_double4 double4;
typedef cl_double3 double3;
typedef cl_double2 double2;

typedef cl_int4 int4;
typedef cl_int3 int3;
typedef cl_int2 int2;

typedef cl_uint4 uint4;
typedef cl_uint3 uint3;
typedef cl_uint2 uint2;

double rsqrt(double x);

#define INT4(_x, _y, _z, _w) (int4)   {.x = (int)(_x),    .y = (int)(_y),    .z = (int)(_z),   .w = (int)(_w)}
#define INT2(_x, _y)         (int2)   {.x = (int)(_x),    .y = (int)(_y)}
#define UINT3(_x, _y, _z)    (uint3)  {.x = (uint)(_x),   .y = (uint)(_y),   .z = (uint)(_z)}
#define UINT2(_x, _y)        (uint2)  {.x = (uint)(_x),   .y = (uint)(_y)}
#define DOUBLE3(_x, _y, _z)  (double3){.x = (double)(_x), .y = (double)(_y), .z = (double)(_z)}

#else

#define UINT2(_x, _y)        (uint2)  ((_x), (_y))
#define UINT3(_x, _y, _z)    (uint3)  ((_x), (_y), (_z))
#define INT2(_x, _y)         (int2)   ((_x), (_y))
#define INT4(_x, _y, _z, _w) (int4)   ((_x), (_y), (_z), (_w))
#define DOUBLE3(_x, _y, _z)  (double3)((_x), (_y), (_z))

typedef unsigned long size_t;

#endif

#endif
