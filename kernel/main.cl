#include <opencl-c-base.h>

#include "bbox.h"
#include "scene.h"

#include "math.c"
#include "float3.c"

kernel void compute_ray(
	__global Sphere *spheres, __global Material *materials,
	__write_only image2d_t output
)
{
	float4 out = (0, 0, 0, 0);
	const int2 pos = {get_global_id(0), get_global_id(1)};
	/*
		w = blue
		x = alpha ?
		y = red
		z = green
	*/

	if(pos.x == spheres[0].center.x && pos.y == spheres[0].center.y)
		out.y = out.z = 0.5f;

	write_imagef(output, pos, out);
}
