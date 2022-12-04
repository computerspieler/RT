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
#define __constant 

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

#define INT4(_x, _y, _z, _w)    (int4)   {.x = (int)(_x),    .y = (int)(_y),    .z = (int)(_z),   .w = (int)(_w)}
#define INT2(_x, _y)            (int2)   {.x = (int)(_x),    .y = (int)(_y)}
#define INT3(_x, _y, _z)        (int3)   {.x = (int)(_x),    .y = (int)(_y),    .z = (int)(_y)}
#define UINT3(_x, _y, _z)       (uint3)  {.x = (uint)(_x),   .y = (uint)(_y),   .z = (uint)(_z)}
#define UINT2(_x, _y)           (uint2)  {.x = (uint)(_x),   .y = (uint)(_y)}
#define DOUBLE2(_x, _y)         (double2){.x = (double)(_x), .y = (double)(_y)}
#define DOUBLE3(_x, _y, _z)     (double3){.x = (double)(_x), .y = (double)(_y), .z = (double)(_z)}
#define FLOAT4(_x, _y, _z, _w)  (float4) {.x = (float)(_x),  .y = (float)(_y),  .z = (float)(_z), .w = (float)(_w)}

#else

#define UINT2(_x, _y)           (uint2)  ((_x), (_y))
#define UINT3(_x, _y, _z)       (uint3)  ((_x), (_y), (_z))
#define INT2(_x, _y)            (int2)   ((_x), (_y))
#define INT3(_x, _y, _z)        (int3)   ((_x), (_y), (_z))
#define INT4(_x, _y, _z, _w)    (int4)   ((_x), (_y), (_z), (_w))
#define DOUBLE2(_x, _y)         (double2)((_x), (_y))
#define DOUBLE3(_x, _y, _z)     (double3)((_x), (_y), (_z))
#define FLOAT4(_x, _y, _z, _w)  (float4) ((_x), (_y), (_z), (_w))

typedef unsigned long size_t;

#endif

#define VEC3(_x, _y, _z)        DOUBLE3(_x, _y, _z)
typedef double3 vec3;

#endif
