#ifndef _CL_HANDLER_H_
#define _CL_HANDLER_H_

#include <stdio.h>
#include <stdbool.h>
#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <CL/opencl.h>

#include "array.h"

typedef struct CLBuffer CLBuffer;
struct CLBuffer
{
    cl_mem cl;
    void *buf;
    size_t len;

    bool input;
    bool output;
    bool image;
};

typedef struct OpenCL_GeneralContext OpenCL_GeneralContext;
struct OpenCL_GeneralContext
{  
    cl_platform_id platform;
    cl_device_id device;

    cl_context context;
    cl_command_queue queue;
};

typedef struct OpenCL_ProgramContext OpenCL_ProgramContext;
struct OpenCL_ProgramContext
{
    cl_program program;
    cl_kernel kern;

    Array buffers;
    Array programs_src;
    Array programs_src_len;
};

#define CHECK_ERR(func, on_err)         \
    if(err != CL_SUCCESS) {             \
        printf(__FILE__ "(%d) " #func ": %d\n", __LINE__, err);    \
        goto on_err;                    \
    }

#define TRY(func, on_err, ...)          \
    err = func(__VA_ARGS__);            \
    CHECK_ERR(func, on_err)

void CL_CALLBACK build_notify (cl_program program, void *user_data);
cl_int opencl_init_general_context(OpenCL_GeneralContext *cl_gen);
cl_int opencl_init_program_context(OpenCL_ProgramContext *cl_prg);
cl_int opencl_add_program_source(OpenCL_ProgramContext *cl_prg, FILE *f);
cl_int opencl_build_program(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, char* kernel_entry, const char *compile_flags);
int opencl_add_input_buffer(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, void *buf, size_t len);
int opencl_add_input_image(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, cl_uchar4 **raw, int width, int height);
int opencl_add_input_output_buffer(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, void *buf, size_t len);
int opencl_add_output_image(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, cl_uchar4 **raw, int width, int height);
void opencl_prerun(OpenCL_ProgramContext *cl_prg);
void opencl_run(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, size_t *origin, size_t *region);
void opencl_postrun(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, size_t *origin, size_t *region);

#endif