#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <SDL/SDL.h>
#include <string.h>

#include "render.h"

const int SCREEN_WIDTH  = 640;
const int SCREEN_HEIGHT = 480;

SDL_Surface* screen;

int initWindow()
{
	if(SDL_Init(SDL_INIT_VIDEO) < 0)
		SDL_FATAL(SDL_Init);

	screen = SDL_SetVideoMode(
			SCREEN_WIDTH,
			SCREEN_HEIGHT,
			32, SDL_HWSURFACE
		);

	if(!screen)
		SDL_FATAL(SDL_SetVideoMode);

	SDL_WM_SetCaption("RayTracing", NULL);

	return 0;
}

int updateWindow(char* output)
{
	SDL_Event ev;
	int running = 1;

    blitBuffer(screen, output, 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

    SDL_Flip(screen);

	while(running)
	{
		if(!SDL_PollEvent(&ev))
			continue;

		switch(ev.type)
		{
			case SDL_QUIT:
				running = 0;
				break;
		}
	}

	return 0;
}

typedef struct {  
    cl_platform_id platform;
    cl_device_id device;

    cl_context context;
    cl_command_queue queue;

    cl_program program;
    cl_kernel kernel;

    cl_mem output;
} OpenCL_FullContext;

// TODO: Better error handling for OpenCL
cl_int opencl_init(OpenCL_FullContext *cl)
{
    cl_int err;

    clGetPlatformIDs(1, &cl->platform, NULL);
    clGetDeviceIDs(cl->platform, CL_DEVICE_TYPE_GPU, 1, &cl->device, NULL);

    cl->context = clCreateContext(NULL, 1, &cl->device, NULL, NULL, &err);
    if(err != CL_SUCCESS)
    {
        printf("clCreateContext: ");
        return err;
    }

    cl->queue = clCreateCommandQueueWithProperties(cl->context, cl->device, NULL, &err);
    if(err != CL_SUCCESS)
        printf("clCreateCommandQueueWithProperties: ");
    
    return err;
}

int read_full_file(FILE *file, char **file_str, size_t *file_str_len)
{
    assert(file);
    assert(file_str);
    assert(file_str_len);

    size_t length_read;

    fseek(file, 0, SEEK_END);
    *file_str_len = ftell(file);
    rewind(file);

    *file_str = (char*) malloc(*file_str_len);
    if(!*file_str)
    {
        perror("malloc");
        return -1;
    }
    
    length_read = fread(*file_str, sizeof(char), *file_str_len, file);
    if(length_read != *file_str_len)
    {
        perror("fread");
        return -1;
    }

    return 0;
}

cl_int opencl_build_program(FILE *program, OpenCL_FullContext *cl, char* kernel_entry)
{
    cl_int err;

    char *program_str;
    size_t program_str_len;

    if(read_full_file(program, &program_str, &program_str_len))
        return -1;
    
    cl->program = clCreateProgramWithSource(cl->context, 1, (const char **) &program_str, &program_str_len, &err);
    if(err != CL_SUCCESS)
    {
        printf("clCreateProgramWithSource: ");
        return err;
    }
    
    err = clBuildProgram(cl->program, 1, &cl->device, NULL, NULL, NULL);
    if(err != CL_SUCCESS)
    {
        printf("clBuildProgram: ");
        return err;
    }

    cl->kernel = clCreateKernel(cl->program, kernel_entry, &err);
    if(err != CL_SUCCESS)
        printf("clCreateKernel: ");
    
    return err;
}

int main(void)
{
    FILE* fi;
    cl_int err;
    char* output = NULL;
    OpenCL_FullContext cl;
    
    printf("Opening the context\n");
    err = opencl_init(&cl);
    if(err != CL_SUCCESS)
        goto exit;

    fi = fopen("kernel.cl", "r");
    if(!fi) {
        perror("fopen");
        goto exit;
    }

    printf("Compiling the program\n");
    err = opencl_build_program(fi, &cl, "compute_ray");
    if(err != CL_SUCCESS)
        goto exit;
    fclose(fi);

    printf("Creating the output buffer\n");
    const cl_image_format fmt = {
        .image_channel_order = CL_RGBA,
        .image_channel_data_type = CL_UNSIGNED_INT8
    };

    cl.output = clCreateImage2D(
        cl.context,
        CL_MEM_WRITE_ONLY,
        &fmt,
        SCREEN_WIDTH, SCREEN_HEIGHT, 0,
        NULL, &err
    );
    if(err != CL_SUCCESS)
        goto exit;

    printf("Here we go\n");
    clSetKernelArg(cl.kernel, 0, sizeof(cl_mem), (void*) &cl.output);

    output = malloc(SCREEN_WIDTH*SCREEN_HEIGHT*4);
    memset(output, 0, SCREEN_WIDTH*SCREEN_HEIGHT*4);

    size_t offset[3] = {0};
    size_t size[3] = {SCREEN_WIDTH, SCREEN_HEIGHT,1};
    err = clEnqueueNDRangeKernel(cl.queue, cl.kernel, 2, offset, size, NULL, 0, NULL, NULL);
    if(err != CL_SUCCESS)
        goto exit;

    size_t origin[3] = {0};
    size_t region[3] = {SCREEN_WIDTH, SCREEN_HEIGHT,1};
    err = clEnqueueReadImage(cl.queue, cl.output, CL_TRUE, origin, region, 0, 0, output, 0, NULL, NULL);
    if(err != CL_SUCCESS)
        goto exit;

    printf("Done !\n");

    if(initWindow())
        goto exit;

    updateWindow(output);

exit:
    printf("Error code : %i\n", err);
    printf("Exit !\n");

    clReleaseMemObject(cl.output);
    clReleaseKernel(cl.kernel);
    clReleaseProgram(cl.program);
    clReleaseCommandQueue(cl.queue);
    clReleaseContext(cl.context);

    if(output)
        free(output);
    return 0;
}