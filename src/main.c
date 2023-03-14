#include <CL/cl.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <SDL/SDL.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

#include "bvhtree.h"
#include "cl_handler.h"
#include "material.h"
#include "scene_loader.h"
#include "scene.h"
#include "camera.h"
#include "transform.h"
#include "host.h"
#include "kernel/main.cl.h"

#define SDL_FATAL(func)						\
{											\
	printf(#func ": %s\n", SDL_GetError());	\
	exit(EXIT_FAILURE);						\
}

#define SCREEN_WIDTH  400
#define SCREEN_HEIGHT 400

#define SURFACE_WIDTH  100
#define SURFACE_HEIGHT 100

long long get_current_time()
{
	struct timeval te;
	gettimeofday(&te, NULL);
	return te.tv_sec * 1000LL + te.tv_usec / 1000;
}

const cl_double speed = 0.125f;
int main(int argc, char* argv[])
{
	bool update;
	FILE *f;
	Scene obj;
	cl_int err;
	SDL_Surface* screen;
	OpenCL_GeneralContext cl_gen;
	OpenCL_ProgramContext cl_prg;
	SDL_Event ev;
	int running = 1;
	int selected_triangle = -1;
	size_t origin[3] = {0};
	size_t region[3] = {SCREEN_WIDTH, SCREEN_HEIGHT, 1};
	void *output;

	HostContext ctx = (HostContext) {
		.sampleCount = 0,
		.camera = (Camera) {
			.near = 0.1f,
			.pos = VEC3(0, 0, 0),
			.rot = VEC3(0, 0, 0),
			.fov = M_PI / 3.f,
			.max_t = 1e2,
			.viewport = (int2){.x = SCREEN_WIDTH, .y = SCREEN_HEIGHT}
		},
		.simpleView = true,
		.max_threshold = 10
	};
	ctx.camera.rotation_transform = transform_combine(
		transform_combine(
			transform_rotate_x(ctx.camera.rot.x),
			transform_rotate_y(ctx.camera.rot.y)
		), transform_rotate_z(ctx.camera.rot.z)
	);

	if(argc != 2) {
		printf("Usage: program [SCENE]\n");
		return -1;
	}

	srand(time(NULL));

	if(SDL_Init(SDL_INIT_VIDEO) < 0)
		SDL_FATAL(SDL_Init);

	if (SDL_EnableKeyRepeat(100, SDL_DEFAULT_REPEAT_INTERVAL))
		SDL_FATAL(SDL_EnableKeyRepeat);

	screen = SDL_SetVideoMode(
		SCREEN_WIDTH,
		SCREEN_HEIGHT,
		32, SDL_HWSURFACE
	);

	if(!screen)
		SDL_FATAL(SDL_SetVideoMode);

	SDL_WM_SetCaption("RayTracing", NULL);
	
	// Open the context
	err = opencl_init_general_context(&cl_gen);
	if(err != CL_SUCCESS)
		goto exit;

	// Load the kernel
	err = opencl_init_program_context(&cl_prg);
	if(err != CL_SUCCESS)
		goto exit;

	err = opencl_add_program_source(&cl_prg, kernel_main_cl_src, sizeof(kernel_main_cl_src));
	if(err != CL_SUCCESS)
		goto exit;

	// Load the scene
	f = fopen(argv[1], "r");
	if(loadFromObj(f, &obj, &ctx.scene_meta, false))
		goto exit;
	fclose(f);

	printObjectInfo(ctx.scene_meta);

	char buffer[1024];

	printf("Compile the program\n");
	sprintf(
		buffer, "-DMAX_TREE_DEPTH=%u",
		BVHTree_depth(ctx.scene_meta.tree.nodes, ctx.scene_meta.tree.root)
	);

	TRY(opencl_build_program, exit,
		&cl_gen, &cl_prg, "compute_pixel", buffer);
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.materials, ctx.scene_meta.materials_count * sizeof(Material));
	
	printf("Load the buffers\n");
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.vertices, ctx.scene_meta.vertices_count * sizeof(vec3));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.vertices_tex, ctx.scene_meta.vertices_tex_count * sizeof(vec2));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.triangles, ctx.scene_meta.triangles_count * sizeof(Triangle));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, ctx.scene_meta.tree.nodes, ctx.scene_meta.tree.nodes_count * sizeof(BVHNode));

	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, &ctx, sizeof(HostContext));

	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, &selected_triangle, sizeof(int));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.lights, ctx.scene_meta.lights_count * sizeof(Light));

	output = malloc(SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(Float));
	bzero(output, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(Float));
	
	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, output, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(Float));

	TRY(opencl_add_output_image, exit,
		&cl_gen, &cl_prg, (cl_uchar4**) &screen->pixels, SCREEN_WIDTH, SCREEN_HEIGHT
	);

	// 
	void *irradiance = malloc(SURFACE_WIDTH * SURFACE_HEIGHT * sizeof(Float));

	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, irradiance, SURFACE_WIDTH * SURFACE_HEIGHT * sizeof(Float));

	printf("Run the program\n");
	opencl_prerun(&cl_prg);
    ctx.lambda = 400;

    printf("Wavelength: %d\n", ctx.lambda);
	long long first_start = get_current_time();
	float mean_time = 0.f;
	
	while(running)
	{
		while(!SDL_PollEvent(&ev) || ev.type == SDL_NOEVENT) {
			SDL_LockSurface(screen);
			for(int i = 0; i < 1; i ++) {
				ctx.sampleCount ++;
				long long start = get_current_time();
				ctx.seed = (uint) start;
				opencl_run(&cl_gen, &cl_prg, origin, region);
				opencl_postrun(&cl_gen, &cl_prg, origin, region);
				mean_time *= ctx.sampleCount - 1;
				mean_time += get_current_time() - start;
				mean_time /= ctx.sampleCount;
				printf(
					"Milliseconds spent: %lld, Mean time per sampling: %f, time since the beginning: %lld, Wavelength: %d\n",
					get_current_time() - start,
					mean_time,
					(get_current_time() - first_start) / 1000,
					ctx.lambda
				);
			}
			SDL_UnlockSurface(screen);
			
			SDL_Flip(screen);
		}
		
		switch(ev.type)
		{
			case SDL_QUIT:
				running = 0;
				break;

			case SDL_KEYDOWN:
				update = ctx.simpleView;
				switch(ev.key.keysym.sym) {
				case SDLK_q:			if(ctx.simpleView) ctx.camera.pos.x -= speed; break;
				case SDLK_d:			if(ctx.simpleView) ctx.camera.pos.x += speed; break;
				case SDLK_SPACE:		if(ctx.simpleView) ctx.camera.pos.y -= speed; break;
				case SDLK_BACKSPACE:	if(ctx.simpleView) ctx.camera.pos.y += speed; break;
				case SDLK_z:			if(ctx.simpleView) ctx.camera.pos.z += speed; break;
				case SDLK_s:			if(ctx.simpleView) ctx.camera.pos.z -= speed; break;

				case SDLK_UP:			if(ctx.simpleView) ctx.camera.rot.x += speed; break;
				case SDLK_DOWN:			if(ctx.simpleView) ctx.camera.rot.x -= speed; break;
				case SDLK_LEFT:			if(ctx.simpleView) ctx.camera.rot.y -= speed; break;
				case SDLK_RIGHT:		if(ctx.simpleView) ctx.camera.rot.y += speed; break;
				case SDLK_PAGEUP:		if(ctx.simpleView) ctx.camera.rot.z -= speed; break;
				case SDLK_PAGEDOWN:		if(ctx.simpleView) ctx.camera.rot.z += speed; break;
				
				case SDLK_p:			ctx.max_threshold += 1;	break;
				case SDLK_m:			ctx.max_threshold -= 1;	break;
				case SDLK_TAB:
					update = true;
					ctx.simpleView ^= 1;
					break;

				default:
					update = false;
				}

				if(!update)
					break;

				mean_time = 0;
				ctx.sampleCount = 0;

				ctx.camera.rotation_transform = transform_combine(
					transform_combine(
						transform_rotate_x(ctx.camera.rot.x),
						transform_rotate_y(ctx.camera.rot.y)
					), transform_rotate_z(ctx.camera.rot.z)
				);
				
				break;
			
			case SDL_MOUSEMOTION:
				ctx.mouse_pos.x = ev.motion.x;
				ctx.mouse_pos.y = ev.motion.y;
				break;

			default:
				break;
		}
	}

exit:
	printf("Error code : %i\n", err);
	printf("Exit !\n");

	return 0;
}
