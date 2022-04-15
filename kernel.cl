
#define RAY_EPSILON 1e-3d

enum ObjectType
{
    SPHERE, PLANE
};

struct Ray
{
	float3 o;				// Origine
	float3 d;				// Direction

	double t;				// Le temps
	double min_t, max_t;	// L'interval dans lequel le temps doit être compris
};

struct Scene {
    enum ObjectType *types;
    int3 *pos;
    void **metadata;
};

kernel void compute_ray(__write_only image2d_t output)
{
    const int2 pos = {get_global_id(0), get_global_id(1)};
    int4 out = 0;
    out.y = pos.x;
    out.z = pos.y;
    write_imagei (output, (int2)(pos.x, pos.y), out);
}
