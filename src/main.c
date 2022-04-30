#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <SDL/SDL.h>
#include <string.h>

#include "util/render.h"
#include "scene.h"

#define CL_COMPILE_FLAGS "-D IN_OPENCL -I include -I kernel -I common"
 
/*
	float3 center;
	float radius;
	int material;
*/

Material sample_mats[] = {
    (Material){
        .color = (float3) {0.5f, 0.5f, 0}
    },
    (Material){
        .color = (float3) {0, 0.5f, 0}
    },
};
Sphere sample_spheres[] = {
    (Sphere){
        .center = (float3) {1, 0, 0},
        .radius = 5,
        .material = 0
    }
};
Scene sample_scene = {
    .materials = sample_mats,
    .materials_count = 1,

    .spheres = sample_spheres,
    .spheres_count = 1
};

typedef struct {  
    cl_platform_id platform;
    cl_device_id device;

    cl_context context;
    cl_command_queue queue;

    cl_program program;
    cl_kernel kernel;
} OpenCL_FullContext;

void CL_CALLBACK build_notify (cl_program program, void *user_data)
{
    char *buffer = NULL;
    size_t buffer_len;
    cl_int err;
    cl_build_status status;
    OpenCL_FullContext *cl = user_data;

    err = clGetProgramBuildInfo(program, cl->device, CL_PROGRAM_BUILD_STATUS, sizeof(status), &status, NULL);
    if(err != CL_SUCCESS)
        printf("clGetProgramBuildInfo: %d\n", err);

    else if(status == CL_BUILD_ERROR) {
        err = clGetProgramBuildInfo(program, cl->device, CL_PROGRAM_BUILD_LOG, 0, NULL, &buffer_len);
        buffer = malloc(buffer_len);
        err = clGetProgramBuildInfo(program, cl->device, CL_PROGRAM_BUILD_LOG, buffer_len, buffer, NULL);
        if(err != CL_SUCCESS)
            printf("clGetProgramBuildInfo(2): %d\n", err);
        else
            printf("Error:\n%.*s\n", (int) buffer_len, buffer);

        free(buffer);
    }
}

// TODO: Better error handling for OpenCL
cl_int opencl_init(OpenCL_FullContext *cl)
{
    cl_int err = CL_SUCCESS;

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
    
    err = clBuildProgram(cl->program, 1, &cl->device, CL_COMPILE_FLAGS, build_notify, cl);
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

int main(int argc, char* argv[])
{
    FILE* fi;
    cl_int err;
    OpenCL_FullContext cl;

    cl_mem output = NULL;
    cl_mem sphere = NULL;
    cl_mem material = NULL;
    char* output_raw = NULL;

    if(argc != 2) {
        printf("Only the filename");
        return -1;
    }
    
    printf("Opening the context\n");
    err = opencl_init(&cl);
    if(err != CL_SUCCESS)
        goto exit;

    fi = fopen(argv[1], "r");
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

    output = clCreateImage2D(
        cl.context,
        CL_MEM_WRITE_ONLY,
        &fmt,
        SCREEN_WIDTH, SCREEN_HEIGHT, 0,
        NULL, &err
    );

    if(err != CL_SUCCESS)
    {
        printf("clCreateImage2D: ");
        goto exit;
    }

    sphere = clCreateBuffer(
        cl.context,
        CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
        sample_scene.spheres_count * sizeof(Sphere),
        sample_scene.spheres,
        &err
    );

    if(err != CL_SUCCESS)
    {
        printf("clCreateBuffer: ");
        goto exit;
    }

    material = clCreateBuffer(
        cl.context,
        CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR,
        sample_scene.materials_count * sizeof(Material),
        sample_scene.materials,
        &err
    );

    if(err != CL_SUCCESS)
    {
        printf("clCreateBuffer: ");
        goto exit;
    }

    printf("Here we go\n");
    clSetKernelArg(cl.kernel, 0, sizeof(cl_mem), (void*) &sphere);
    clSetKernelArg(cl.kernel, 1, sizeof(cl_mem), (void*) &material);
    clSetKernelArg(cl.kernel, 2, sizeof(cl_mem), (void*) &output);

    output_raw = malloc(SCREEN_WIDTH*SCREEN_HEIGHT*4);
    memset(output_raw, 0, SCREEN_WIDTH*SCREEN_HEIGHT*4);

    // Run the program
    size_t offset[3] = {0};
    size_t size[3] = {SCREEN_WIDTH, SCREEN_HEIGHT,1};
    err = clEnqueueNDRangeKernel(cl.queue, cl.kernel, 2, offset, size, NULL, 0, NULL, NULL);
    if(err != CL_SUCCESS)
    {
        printf("clEnqueueNDRangeKernel: ");
        goto exit;
    }

    err = clEnqueueWriteBuffer(
        cl.queue, sphere, CL_FALSE, 0,
        sample_scene.spheres_count * sizeof(Sphere),
        sample_scene.spheres,
        0,
        NULL, NULL
    );
    if(err != CL_SUCCESS)
    {
        printf("clEnqueueWriteBuffer: ");
        goto exit;
    }

    err = clEnqueueWriteBuffer(
        cl.queue, material, CL_FALSE, 0,
        sample_scene.materials_count * sizeof(Material),
        sample_scene.materials,
        0,
        NULL, NULL
    );
    if(err != CL_SUCCESS)
    {
        printf("clEnqueueWriteBuffer: ");
        goto exit;
    }

    size_t origin[3] = {0};
    size_t region[3] = {SCREEN_WIDTH, SCREEN_HEIGHT,1};
    err = clEnqueueReadImage(cl.queue, output, CL_TRUE, origin, region, 0, 0, output_raw, 0, NULL, NULL);
    if(err != CL_SUCCESS)
    {
        printf("clEnqueueReadImage: ");
        goto exit;
    }

    printf("Done !\n");

    if(initWindow())
        goto exit;

    updateWindow(output_raw);

exit:
    printf("Error code : %i\n", err);
    printf("Exit !\n");

    clReleaseMemObject(output);
    clReleaseKernel(cl.kernel);
    clReleaseProgram(cl.program);
    clReleaseCommandQueue(cl.queue);
    clReleaseContext(cl.context);

    if(output_raw)
        free(output_raw);
    return 0;
}