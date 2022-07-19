#include <CL/cl.h>
#include <linux/limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <SDL/SDL.h>
#include <string.h>
#include <limits.h>

#include "bvhtree.h"
#include "cl_handler.h"
#include "material.h"
#include "obj_loader.h"
#include "object.h"
#include "render.h"
#include "camera.h"
#include "scene.h"

const char* CL_COMPILE_FLAGS = "-D IN_OPENCL -I include -I kernel -I common";

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

int main(int argc, char *argv[])
{
    FILE *f;
    Object obj;
    ObjectMetadata obj_meta;
    cl_int err;
    OpenCL_GeneralContext cl_gen;
    OpenCL_ProgramContext cl_prg;
    
    f = fopen("scene/test2.obj", "r");
    if(loadFromObj(f, &obj, &obj_meta))
        return -1;

    printObjectInfo(obj_meta);

    Camera camera = {
        .near = 0.1f,
        .pos = DOUBLE3(0, 0, -10),
        .fov = M_PI_4,
        .max_t = 1e10,
        .viewport = (int2){.x = SCREEN_WIDTH, .y = SCREEN_HEIGHT}
    };
    
    // Open the context
    err = opencl_init_general_context(&cl_gen);
    if(err != CL_SUCCESS)
        goto exit;

    // Load the kernel
    err = opencl_init_program_context(&cl_prg);
    if(err != CL_SUCCESS)
        goto exit;

    for(int i = 1; i < argc; i++) {
        err = add_folder(&cl_prg, argv[i]);
        // If this isn't a folder
        if(err == -1) {
            f = fopen(argv[i], "r");
            if(!f) {
                perror("fopen");
                goto exit;
            }

            opencl_add_program_source(&cl_prg, f);
            fclose(f);
        } else if(err != 0)
            goto exit;
    }

    char buffer[1024];
    sprintf(buffer, "%s -DMAX_TREE_DEPTH=%d", CL_COMPILE_FLAGS, BVHTree_depth(obj_meta.tree.nodes, obj_meta.tree.root));

    TRY(opencl_build_program, exit,
        &cl_gen, &cl_prg, "compute_ray", buffer);
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, &camera, sizeof(Camera));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj.materials, obj_meta.materials_count * sizeof(Material));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj.vertices, obj_meta.vertices_count * sizeof(double3));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj.vertices_normal, obj_meta.vertices_normal_count * sizeof(double3));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj.vertices_tex, obj_meta.vertices_tex_count * sizeof(double2));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj.triangles, obj_meta.triangles_count * sizeof(Triangle));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, &obj_meta, sizeof(ObjectMetadata));
    
    TRY(opencl_add_input_buffer, exit,
        &cl_gen, &cl_prg, obj_meta.tree.nodes, obj_meta.tree.nodes_count * sizeof(BVHNode));

    if(initWindow())
        goto exit;

    updateWindow(&cl_gen, &cl_prg, &camera);

exit:
    printf("Error code : %i\n", err);
    printf("Exit !\n");

    return 0;
}
