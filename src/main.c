#include <CL/cl.h>
#include <linux/limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <SDL/SDL.h>
#include <string.h>
#include <limits.h>
#include <bits/types/struct_timeval.h>
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

#define SDL_FATAL(func)						\
{											\
	printf(#func ": %s\n", SDL_GetError());	\
	exit(EXIT_FAILURE);						\
}

#define SCREEN_WIDTH  320
#define SCREEN_HEIGHT 320

#define SURFACE_WIDTH  100
#define SURFACE_HEIGHT 100

const char* CL_COMPILE_FLAGS =
	" -D IN_OPENCL "
	" -I include "
	" -I kernel "
	" -I common "
	" -D PRINT_DEBUG "
	;

cl_int add_folder(OpenCL_ProgramContext *cl, const char *path)
{
	DIR *d;
	FILE *f;
	struct dirent *dir;
	char fullpath[PATH_MAX];

	d = opendir(path);
	if(!d) {
		perror("opendir");
		return -1;
	}

	while((dir = readdir(d)) != NULL) {
		if(dir->d_type != DT_REG)
			continue;

		sprintf(fullpath, "%s/%s", path, dir->d_name);
		printf("Load %s\n", fullpath);

		f = fopen(fullpath, "r");
		if(!f) {
			perror("fopen");
			closedir(d);
			return -2;
		}

		opencl_add_program_source(cl, f);
		fclose(f);
	}

	closedir(d);

	return CL_SUCCESS;
}

long long get_current_time()
{
	struct timeval te;
	gettimeofday(&te, NULL);
	return te.tv_sec * 1000LL + te.tv_usec / 1000;
}

const cl_double speed = 0.125f;
int main(void)
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
	size_t origin[3] = {0};
	size_t region[3] = {SCREEN_WIDTH, SCREEN_HEIGHT, 1};
	void *b;

	HostContext ctx = (HostContext) {
		.camera = (Camera) {
			.near = 0.1f,
			.pos = VEC3(0, 0, 0),
			.rot = VEC3(0, 0, 0),
			.fov = M_PI_4,
			.max_t = 1e2,
			.viewport = (int2){.x = SCREEN_WIDTH, .y = SCREEN_HEIGHT}
		},
		.simpleView = true,
		.selected_triangle = -1
	};
	ctx.camera.rotation_transform = transform_combine(
		transform_combine(
			transform_rotate_x(ctx.camera.rot.x),
			transform_rotate_y(ctx.camera.rot.y)
		), transform_rotate_z(ctx.camera.rot.z)
	);

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

	err = add_folder(&cl_prg, "common");
	if(err != CL_SUCCESS)
		goto exit;

	err = add_folder(&cl_prg, "kernel");
	if(err != CL_SUCCESS)
		goto exit;

	// Load the scene
	f = fopen("scene/Room.obj", "r");
	if(loadFromObj(f, &obj, &ctx.scene_meta, false))
		goto exit;
	fclose(f);

	printObjectInfo(ctx.scene_meta);

	char buffer[1024];

	printf("Compile the program\n");
	sprintf(
		buffer, "%s -DMAX_TREE_DEPTH=%u",
		CL_COMPILE_FLAGS,
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
		&cl_gen, &cl_prg, obj.vertices_tex, ctx.scene_meta.vertices_tex_count * sizeof(double2));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.triangles, ctx.scene_meta.triangles_count * sizeof(Triangle));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, ctx.scene_meta.tree.nodes, ctx.scene_meta.tree.nodes_count * sizeof(BVHNode));

	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, &ctx, sizeof(HostContext));

	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.map, ctx.scene_meta.map_size * sizeof(unsigned char));
	
	b = malloc(SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(double));
	bzero(b, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(double));
	
	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, b, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(double));

	TRY(opencl_add_output_image, exit,
		&cl_gen, &cl_prg, (cl_uchar4**) &screen->pixels, SCREEN_WIDTH, SCREEN_HEIGHT
	);

	// 
	void *irradiance = malloc(SURFACE_WIDTH * SURFACE_HEIGHT * sizeof(double));

	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, irradiance, SURFACE_WIDTH * SURFACE_HEIGHT * sizeof(double));

	printf("Run the program\n");
	ctx.sampleCount = 0;
	opencl_prerun(&cl_prg);
    ctx.lambda = 400;

    printf("Wavelength: %d\n", ctx.lambda);
	long long first_start = get_current_time();
	while(running)
	{
		while(!SDL_PollEvent(&ev) || ev.type == SDL_NOEVENT) {
			SDL_LockSurface(screen);
			for(int i = 0; i < 1; i ++) {
				long long start = get_current_time();
				ctx.seed = (unsigned int) get_current_time() * get_current_time();
				ctx.sampleCount++;
				opencl_run(&cl_gen, &cl_prg, origin, region);
				opencl_postrun(&cl_gen, &cl_prg, origin, region);
				printf("Milliseconds spent: %lld, sampling count: %d, time since the beginning: %lld, Wavelength: %d\n",
					get_current_time() - start, ctx.sampleCount, (get_current_time() - first_start) / 1000,
					ctx.lambda);
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
				update = true;
				switch(ev.key.keysym.sym) {
				case SDLK_q:			ctx.camera.pos.x -= speed; break;
				case SDLK_d:			ctx.camera.pos.x += speed; break;
				case SDLK_SPACE:		ctx.camera.pos.y -= speed; break;
				case SDLK_BACKSPACE:	ctx.camera.pos.y += speed; break;
				case SDLK_z:			ctx.camera.pos.z += speed; break;
				case SDLK_s:			ctx.camera.pos.z -= speed; break;

				case SDLK_UP:			ctx.camera.rot.x += speed; break;
				case SDLK_DOWN:			ctx.camera.rot.x -= speed; break;
				case SDLK_LEFT:			ctx.camera.rot.y -= speed; break;
				case SDLK_RIGHT:		ctx.camera.rot.y += speed; break;
				case SDLK_PAGEUP:		ctx.camera.rot.z -= speed; break;
				case SDLK_PAGEDOWN:		ctx.camera.rot.z += speed; break;
				
				case SDLK_p:			ctx.lambda ++;			   break;
				case SDLK_m:			ctx.lambda --;			   break;
				case SDLK_TAB:			ctx.simpleView ^= 1;	   break;

				default:
					update = false;
				}

				if(!update)
					break;

				ctx.camera.rotation_transform = transform_combine(
					transform_combine(
						transform_rotate_x(ctx.camera.rot.x),
						transform_rotate_y(ctx.camera.rot.y)
					), transform_rotate_z(ctx.camera.rot.z)
				);
				
				bzero(b, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(double) * 4);
				ctx.sampleCount = 0;
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
