#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <assert.h>
#include <stddef.h>

#include "array.h"
#include "cl_handler.h"

int read_full_file(FILE *file, char **file_str, size_t *file_str_len);

void CL_CALLBACK build_notify (cl_program program, void *user_data)
{
    char *buffer = NULL;
    size_t buffer_len;
    cl_int err;
    cl_build_status status;
    OpenCL_GeneralContext *cl = user_data;

    err = clGetProgramBuildInfo(program, cl->device, CL_PROGRAM_BUILD_STATUS, sizeof(status), &status, NULL);
    CHECK_ERR(clGetProgramBuildInfo, end);

    if(status == CL_BUILD_ERROR) {
        err = clGetProgramBuildInfo(program, cl->device, CL_PROGRAM_BUILD_LOG, 0, NULL, &buffer_len);
        buffer = malloc(buffer_len);
        TRY(clGetProgramBuildInfo, endl,
            program, cl->device, CL_PROGRAM_BUILD_LOG, buffer_len, buffer, NULL);

        printf("=== Build error ===\n%.*s\n", (int) buffer_len, buffer);
endl:
        free(buffer);
    }
end:
    return;
}

cl_int opencl_init_general_context(OpenCL_GeneralContext *cl_gen)
{
    cl_int err = CL_SUCCESS;

    err = clGetPlatformIDs(1, &cl_gen->platform, NULL);
    CHECK_ERR(clGetPlatformIDs, end);
    err = clGetDeviceIDs(cl_gen->platform, CL_DEVICE_TYPE_GPU, 1, &cl_gen->device, NULL);
    CHECK_ERR(clGetDeviceIDs, end);
    cl_gen->context = clCreateContext(NULL, 1, &cl_gen->device, NULL, NULL, &err);
    CHECK_ERR(clCreateContext, end);
    cl_gen->queue = clCreateCommandQueueWithProperties(cl_gen->context, cl_gen->device, NULL, &err);
    CHECK_ERR(clCreateCommandQueueWithProperties, end);
end:
    return err;
}
cl_int opencl_init_program_context(OpenCL_ProgramContext *cl_prg)
{
    cl_prg->buffers = array_create(sizeof(CLBuffer));
    cl_prg->programs_src = array_create(sizeof(char*));
    cl_prg->programs_src_len = array_create(sizeof(size_t));

    return CL_SUCCESS;
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

cl_int opencl_add_program_source(OpenCL_ProgramContext *cl_prg, FILE *f)
{
    char *program_str;
    size_t program_str_len;

    if(read_full_file(f, 
        &program_str, &program_str_len))
        return -1;
    
    array_push(&cl_prg->programs_src, &program_str);
    array_push(&cl_prg->programs_src_len, &program_str_len);

    return CL_SUCCESS;
}

cl_int opencl_build_program(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, char* kernel_entry, const char *compile_flags)
{
    cl_int err;
    
    cl_prg->program = clCreateProgramWithSource(cl_gen->context, 
        array_size(&cl_prg->programs_src), 
        (const char **) cl_prg->programs_src.array, 
        cl_prg->programs_src_len.array,
        &err
    );
    CHECK_ERR(clCreateProgramWithSource, end);

    TRY(clBuildProgram, end,
        cl_prg->program, 1, &cl_gen->device, compile_flags, build_notify, cl_gen
    );

    cl_prg->kern = clCreateKernel(cl_prg->program, kernel_entry, &err);
    CHECK_ERR(clCreateKernel, end);

end:  
    return err;
}

int opencl_add_input_buffer(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, void *buf, size_t len)
{
    cl_int err;
    CLBuffer buffer;
    
    buffer = (CLBuffer) {
        .cl = clCreateBuffer(
            cl_gen->context,
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            len, buf,
            &err
        ),
        .buf = buf,
        .len = len,
        .image = false,
        .input = true,
        .output = false
    };
    CHECK_ERR(clCreateBuffer, exit);    
    array_push(&cl_prg->buffers, &buffer);

exit:
    return err;
}

int opencl_add_input_output_buffer(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, void *buf, size_t len)
{
    cl_int err;
    CLBuffer buffer;
    
    buffer = (CLBuffer) {
        .cl = clCreateBuffer(
            cl_gen->context,
            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
            len, buf,
            &err
        ),
        .buf = buf,
        .len = len,
        .image = false,
        .input = false,
        .output = true
    };
    CHECK_ERR(clCreateBuffer, exit);    
    array_push(&cl_prg->buffers, &buffer);

exit:
    return err;
}

int opencl_add_output_image(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, cl_uchar4 **raw, int width, int height)
{
    cl_int err;

    if(!*raw)
        *raw = malloc(width * height * sizeof(cl_uchar4));

    cl_image_desc output_buffer_desc = { 0 };
    output_buffer_desc.image_width  = width;
    output_buffer_desc.image_height = height;
    output_buffer_desc.image_type   = CL_MEM_OBJECT_IMAGE2D;

    const cl_image_format fmt = {
        .image_channel_order = CL_RGBA,
        .image_channel_data_type = CL_UNSIGNED_INT8
    };
    
    CLBuffer buffer = {
        .cl = clCreateImage(
            cl_gen->context, CL_MEM_WRITE_ONLY,
            &fmt, &output_buffer_desc,
            NULL, &err
        ),
        .image = true,
        .input = false,
        .output = true,
        .buf = *raw,
        .len = width * height * sizeof(cl_uchar4)
    };
    CHECK_ERR(clCreateImage, exit);
    array_push(&cl_prg->buffers, &buffer);

exit:
    return err;   
}

int opencl_add_input_image(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, cl_uchar4 **raw, int width, int height)
{
    cl_int err;

    if(!*raw)
        *raw = malloc(width * height * sizeof(cl_uchar4));

    cl_image_desc output_buffer_desc = { 0 };
    output_buffer_desc.image_width  = width;
    output_buffer_desc.image_height = height;
    output_buffer_desc.image_type   = CL_MEM_OBJECT_IMAGE2D;

    const cl_image_format fmt = {
        .image_channel_order = CL_RGBA,
        .image_channel_data_type = CL_UNSIGNED_INT8
    };
    
    CLBuffer buffer = {
        .cl = clCreateImage(
            cl_gen->context, CL_MEM_READ_ONLY,
            &fmt, &output_buffer_desc,
            NULL, &err
        ),
        .image = true,
        .input = true,
        .output = false,
        .buf = *raw,
        .len = width * height * sizeof(cl_uchar4)
    };
    CHECK_ERR(clCreateImage, exit);
    array_push(&cl_prg->buffers, &buffer);

exit:
    return err;   
}

void opencl_prerun(OpenCL_ProgramContext *cl_prg)
{
    size_t i = 0;
    for(CLBuffer* it = cl_prg->buffers.array; i < array_size(&cl_prg->buffers); it ++, i++) {
        clSetKernelArg(cl_prg->kern, i, sizeof(cl_mem), &it->cl);
    }
}

void opencl_run(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, size_t *origin, size_t *region)
{
    cl_int err;
    size_t local[3] = {8, 8, 1};

    TRY(clEnqueueNDRangeKernel, exit,
        cl_gen->queue, cl_prg->kern, 2, origin, region, local, 0, NULL, NULL);

exit:
    return;
}

void opencl_postrun(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, size_t *origin, size_t *region)
{
    size_t i;
    cl_int err;

    i = 0;
    for(CLBuffer* it = cl_prg->buffers.array; i < array_size(&cl_prg->buffers); it ++, i++) {
        if(it->image) {
            if(it->input)
                TRY(clEnqueueWriteImage, exit,
                    cl_gen->queue, it->cl, CL_TRUE,
                    origin, region,
                    0, 0,
                    it->buf, 0,
                    NULL, NULL
                );
            if(it->output)
                TRY(clEnqueueReadImage, exit,
                    cl_gen->queue, it->cl, CL_TRUE,
                    origin, region,
                    0, 0,
                    it->buf, 0,
                    NULL, NULL
                );
        } else {
            if(it->input) {
                TRY(clEnqueueWriteBuffer, exit,
                    cl_gen->queue, it->cl, CL_FALSE, 0,
                    it->len, it->buf,
                    0, NULL, NULL
                );
            }
            
            if(it->output) {
                TRY(clEnqueueReadBuffer, exit,
                    cl_gen->queue, it->cl, CL_TRUE,
                    0, it->len,
                    it->buf, 0,
                    NULL, NULL
                );
            }
        }
    }

exit:
    return;
}
