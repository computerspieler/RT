#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "bbox.h"
#include "double3.h"
#include "typedef.h"
#include "material.h"
#include "obj_loader.h"
#include "array.h"
#include "object.h"
#include "tree.h"

int loadFromMtl(FILE *f, Node *materials_tree, Array *material_name, Array *materials)
{
    Material mat;
    bool initialized = false;

    int res;
    char *name;
    uint str_length = 0;
    char buffer[1024];

    while(true) {
        res = fscanf(f, "%s", buffer);
        if(res == EOF)
            break;
        
        if(strcmp(buffer, "newmtl") == 0) {
            if(initialized)
                array_push(materials, &mat);
            
            fscanf(f, "%s\n", buffer);
            str_length = strlen(buffer);
            name = calloc(str_length + 1, sizeof(char));
            memcpy(name, buffer, str_length * sizeof(char));
            name[str_length] = 0;
            array_push(material_name, &name);
            add_word(materials_tree, name, material_name->logical_size);
            mat = (Material) {0};
            if(!initialized)
                initialized = true;
        } else if(strcmp(buffer, "Ns") == 0)
            fscanf(f, "%lf\n", &mat.specular);
        else if(strcmp(buffer, "Ni") == 0)
            fscanf(f, "%lf\n", &mat.density);
        else if(strcmp(buffer, "Ka") == 0)
            fscanf(f, "%lf %lf %lf\n",
                &mat.ambient_color.x, &mat.ambient_color.y, &mat.ambient_color.z);
        else if(strcmp(buffer, "Kd") == 0)
            fscanf(f, "%lf %lf %lf\n",
                &mat.diffuse_color.x, &mat.diffuse_color.y, &mat.diffuse_color.z);
        else if(strcmp(buffer, "Ks") == 0)
            fscanf(f, "%lf %lf %lf\n",
                &mat.specular_color.x, &mat.specular_color.y, &mat.specular_color.z);
        else if(strcmp(buffer, "Ke") == 0)
            fscanf(f, "%lf %lf %lf\n",
                &mat.emissive_color.x, &mat.emissive_color.y, &mat.emissive_color.z);
    }

    if(initialized)
        array_push(materials, &mat);

    return 0;
}

int loadFromObj(FILE *f, Object *obj, ObjectMetadata *metadata)
{
    FILE *fm;
    double3 vertex;
    double2 vertex_uv;
    Triangle triangle;
    BBox3 bounds;

    int res;
    char buffer[1024];
    char *name;
    uint current_material = -1;
    uint str_length = 0;
    Array vertices, tex, triangles, groups, normals, material_name, materials;
    Node *materials_tree, *n;

    assert(obj);
    assert(metadata);

    vertices = array_create(sizeof(double3));
    normals = array_create(sizeof(double3));
    tex = array_create(sizeof(double2));
    groups = array_create(sizeof(char*));
    triangles = array_create(sizeof(Triangle));
    materials = array_create(sizeof(Material));
    materials_tree = malloc(sizeof(Node));
    bzero(materials_tree, sizeof(Node));
    material_name = array_create(sizeof(char*));

    if(!f) {
        perror("File");
        return -1;
    }

    while(true) {
        res = fscanf(f, "%s", buffer);
        if(res == EOF)
            break;

        if(strcmp(buffer, "v") == 0) {
            fscanf(f, "%lf %lf %lf\n", &vertex.x, &vertex.y, &vertex.z);

            if(array_size(&vertices) == 0) {
                bounds = (BBox3) {
                    .max = vertex,
                    .min = vertex
                };
            } else
                bounds = bbox_p_union(bounds, vertex);

            array_push(&vertices, &vertex);
        } else if(strcmp(buffer, "vn") == 0) {
            fscanf(f, "%lf %lf %lf\n", &vertex.x, &vertex.y, &vertex.z);
            array_push(&normals, &vertex);
        } else if(strcmp(buffer, "vt") == 0) {
            fscanf(f, "%lf %lf\n", &vertex_uv.x, &vertex_uv.y);
            array_push(&tex, &vertex_uv);
        } else if(strcmp(buffer, "f") == 0) {
            res = fscanf(f, "%d/%d/%d %d/%d/%d %d/%d/%d\n",
                &triangle.vertices[0], &triangle.uv[0], &triangle.normals[0],
                &triangle.vertices[1], &triangle.uv[1], &triangle.normals[1],
                &triangle.vertices[2], &triangle.uv[2], &triangle.normals[2]);

            triangle.material = current_material;
            triangle.group = array_size(&groups) - 1;
            
            triangle.vertices[0] --;
            triangle.vertices[1] --;
            triangle.vertices[2] --;
            
            triangle.normals[0] --;
            triangle.normals[1] --;
            triangle.normals[2] --;
            
            triangle.uv[0] --;
            triangle.uv[1] --;
            triangle.uv[2] --;

            if(res != 9) {
                printf("Invalid face (%d)\n", res);
                goto err;
            }

            array_push(&triangles, &triangle);
        } else if(strcmp(buffer, "g") == 0) {
            fscanf(f, "%s\n", buffer);
            str_length = strlen(buffer);
            name = calloc(str_length + 1, sizeof(char));
            memcpy(name, buffer, (str_length + 1) * sizeof(char));
            name[str_length] = 0;
            array_push(&groups, &name);
        } else if(strcmp(buffer, "usemtl") == 0) {
            fscanf(f, "%s\n", buffer);
            n = get_node(materials_tree, buffer);

            if(!n) {
                printf("Invalid material %s\n", buffer);
                goto err;
            }
            
            current_material = n->value;
        } else if(strcmp(buffer, "mtllib") == 0) {
            fscanf(f, "%s\n", buffer);
            fm = fopen(buffer, "r");
            if(!fm) {
                perror("fopen (mtllib)");
                goto err;
            }

            loadFromMtl(fm, materials_tree, &material_name, &materials);
            fclose(fm);
        }
    }

    free_tree(materials_tree);

    metadata->groups_count = array_size(&groups);
    metadata->groups_name = groups.array;
    metadata->materials_count = array_size(&materials);
    metadata->triangles_count = array_size(&triangles);
    metadata->vertices_tex_count = array_size(&tex);
    metadata->vertices_count = array_size(&vertices);
    metadata->vertices_normal_count = array_size(&normals);
    metadata->bounds = bounds;

    obj->triangles = triangles.array;
    obj->vertices = vertices.array;
    obj->vertices_tex = tex.array;
    obj->materials = materials.array;
    obj->vertices_normal = normals.array;

    printf("Object file sucessfully loaded\n");

    return 0;

err:
    free_tree(materials_tree);
    array_free(&vertices);
    array_free(&tex);
    array_free(&triangles);
    array_free(&groups);

    return -1;
}

int mergeObjects(Object *obj1, ObjectMetadata *obj1_meta, Object *obj2, ObjectMetadata *obj2_meta, Object *output, ObjectMetadata *output_meta)
{
    assert(obj1);
    assert(obj2);
    assert(obj1_meta);
    assert(obj2_meta);
    assert(output);
    assert(output_meta);

    output_meta->groups_count = obj1_meta->groups_count + obj2_meta->groups_count;
    output_meta->vertices_count = obj1_meta->vertices_count + obj2_meta->vertices_count;
    output_meta->materials_count = obj1_meta->materials_count + obj2_meta->materials_count;
    output_meta->triangles_count = obj1_meta->triangles_count + obj2_meta->triangles_count;
    output_meta->vertices_tex_count = obj1_meta->vertices_tex_count + obj2_meta->vertices_tex_count;
    output_meta->vertices_normal_count = obj1_meta->vertices_normal_count + obj2_meta->vertices_normal_count;
    output_meta->bounds = bbox_b_union(obj1_meta->bounds, obj2_meta->bounds);

    output->materials = malloc(sizeof(Material) * output_meta->materials_count);
    output->triangles = malloc(sizeof(Triangle) * output_meta->triangles_count);
    output->vertices = malloc(sizeof(double3) * output_meta->vertices_count);
    output->vertices_tex = malloc(sizeof(double2) * output_meta->vertices_tex_count);
    output->vertices_normal = malloc(sizeof(double3) * output_meta->vertices_normal_count);
    output_meta->groups_name = malloc(sizeof(char*) * output_meta->groups_count);

    memcpy(output_meta->groups_name, obj1_meta->groups_name, sizeof(char*) * obj1_meta->groups_count);
    memcpy(output_meta->groups_name + obj1_meta->groups_count, obj2_meta->groups_name, sizeof(char*) * obj2_meta->groups_count);

    memcpy(output->materials, obj1->materials, sizeof(Material) * obj1_meta->materials_count);
    memcpy(output->materials + obj1_meta->materials_count, obj2->materials, sizeof(Material) * obj2_meta->materials_count);

    memcpy(output->vertices, obj1->vertices, sizeof(double3) * obj1_meta->vertices_count);
    memcpy(output->vertices + obj1_meta->vertices_count, obj2->vertices, sizeof(double3) * obj2_meta->vertices_count);

    memcpy(output->vertices_tex, obj1->vertices_tex, sizeof(double2) * obj1_meta->vertices_tex_count);
    memcpy(output->vertices_tex + obj1_meta->vertices_tex_count, obj2->vertices_tex, sizeof(double2) * obj2_meta->vertices_tex_count);

    memcpy(output->triangles, obj1->triangles, sizeof(Triangle) * obj1_meta->triangles_count);
    for(size_t i = 0; i < obj2_meta->triangles_count; i++)
        output->triangles[i + obj1_meta->triangles_count] = (Triangle) {
            .group = obj2->triangles[i].group + obj1_meta->groups_count,
            .material = obj2->triangles[i].material + obj1_meta->materials_count,
            .normals = {
                obj2->triangles[i].normals[0] + obj1_meta->vertices_normal_count,
                obj2->triangles[i].normals[1] + obj1_meta->vertices_normal_count,
                obj2->triangles[i].normals[2] + obj1_meta->vertices_normal_count
            },
            .vertices = {
                obj2->triangles[i].vertices[0] + obj1_meta->vertices_count,
                obj2->triangles[i].vertices[1] + obj1_meta->vertices_count,
                obj2->triangles[i].vertices[2] + obj1_meta->vertices_count
            },
            .uv = {
                obj2->triangles[i].uv[0] + obj1_meta->vertices_tex_count,
                obj2->triangles[i].uv[1] + obj1_meta->vertices_tex_count,
                obj2->triangles[i].uv[2] + obj1_meta->vertices_tex_count
            }
        };

    return 0;
}