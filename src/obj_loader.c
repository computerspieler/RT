#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "bbox.h"
#include "bvhtree.h"
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

    obj->triangles = triangles.array;
    obj->vertices = vertices.array;
    obj->vertices_tex = tex.array;
    obj->materials = materials.array;
    obj->vertices_normal = normals.array;

    for(size_t i = 0; i < metadata->triangles_count; i++) {
        double3 points[3] = {
            obj->vertices[obj->triangles[i].vertices[0]],
            obj->vertices[obj->triangles[i].vertices[1]],
            obj->vertices[obj->triangles[i].vertices[2]]
        };
        obj->triangles[i].surfaceNormal = double3_cross(double3_diff(points[1], points[0]), double3_diff(points[2], points[0]));
        obj->triangles[i].surfaceNormal = double3_normalize(obj->triangles[i].surfaceNormal);
        obj->triangles[i].det = -double3_dot(obj->triangles[i].surfaceNormal, points[0]);
    }

    buildBVHTreeFromObject(*obj, metadata);

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

/*
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
            },
            .surfaceNormal = obj2->triangles[i].surfaceNormal,
            .det = obj2->triangles[i].det
        };

    return 0;
}
*/

#define BUCKET_COUNT 10

typedef struct PrimitiveInfo PrimitiveInfo;
struct PrimitiveInfo
{
    BBox3 bounds;
    double3 center;
};

int buildBVHTreeSAH(Object obj, ObjectMetadata *obj_meta, Array *nodes, PrimitiveInfo *primitives_info, Triangle *tris, int start, int end)
{
    int dim;
    BVHNode output;
    BBox3 bounds, centroidBounds;

    if(end <= start)
        return -1;

    // If it's only one triangle
    if((end - start) == 1) {
        output = (BVHNode) {
            .bounds = primitives_info[start].bounds,
            .triangle_start = start,
            .triangle_end = end,
            .sons.x = -1,
            .sons.y = -1
        };
        array_push(nodes, &output);
        return array_size(nodes) - 1;
    }

    bounds = primitives_info[start].bounds;
    for(int i = start+1; i < end; i++)
        bounds = bbox_b_union(bounds, primitives_info[i].bounds);
    
    centroidBounds = (BBox3) {
        .min = primitives_info[start].center,
        .max = primitives_info[start].center
    };
    for(int i = start+1; i < end; i++)
        centroidBounds = bbox_p_union(bounds, primitives_info[i].center);

    dim = bbox_maximum_extent(centroidBounds);

    // If it's a bunch of BBox on the same spot, there's no possible splitting
    if(centroidBounds.max.s[dim] == centroidBounds.min.s[dim]) {
        output = (BVHNode) {
            .bounds = bounds,
            .triangle_start = start,
            .triangle_end = end,
            .sons.x = -1,
            .sons.y = -1
        };
        array_push(nodes, &output);
        return array_size(nodes) - 1;
    }

    // Compute the different buckets
    typedef struct BucketInfo BucketInfo;
    struct BucketInfo {
        int count;
        BBox3 bounds;
        bool init;
    };
    BucketInfo buckets[BUCKET_COUNT];
    for(int i = 0; i < BUCKET_COUNT; i++) 
        buckets[i].init = false;

    for(int i = start; i < end; i++) {
        double3 offset = bbox_p_offset(bounds, primitives_info[i].center);
        int b = BUCKET_COUNT * offset.s[dim];
        if(b == BUCKET_COUNT)
            b = BUCKET_COUNT - 1;
        
        if(!buckets[b].init) {
            buckets[b].bounds = primitives_info[i].bounds;
            buckets[b].count = 0;
            buckets[b].init = true;
        } else
            buckets[b].bounds = bbox_b_union(buckets[b].bounds, primitives_info[i].bounds);

        buckets[b].count ++;
    }

    // Compute the cost
    double cost[BUCKET_COUNT - 1];
    int count0[BUCKET_COUNT - 1];
    for(int i = 0; i < BUCKET_COUNT - 1; i ++) {
        BBox3 b0, b1;
        int count1;
        bool init0, init1;

        init0 = init1 = false;
        for(int j = 0; j <= i; j++) {
            if(!buckets[j].init)
                continue;
            
            if(!init0) {
                b0 = buckets[j].bounds;
                count0[i] = buckets[j].count;
                init0 = true;
            } else {
                b0 = bbox_b_union(b0, buckets[j].bounds);
                count0[i] += buckets[j].count;
            }
        }
        for(int j = i + 1; j < BUCKET_COUNT; j++) {
            if(!buckets[j].init)
                continue;
            
            if(!init1) {
                b1 = buckets[j].bounds;
                count1 = buckets[j].count;
                init1 = true;
            } else {
                b1 = bbox_b_union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }
        }

        if(init0 && init1)
            cost[i] = 0.125 + (double) (count0[i] * bbox_surface_area(b0) + count1 * bbox_surface_area(b1)) / bbox_surface_area(bounds);
        else
            cost[i] = +INFINITY;
    }

    // Get the smallest cost
    int minCostIndex = 0;
    for(int i = 1; i < BUCKET_COUNT - 1; i ++)
        if(cost[i] < cost[minCostIndex])
            minCostIndex = i;

    if(cost[minCostIndex] == +INFINITY) {
        output = (BVHNode) {
            .bounds = bounds,
            .triangle_start = start,
            .triangle_end = end,
            .sons.x = -1,
            .sons.y = -1
        };
        array_push(nodes, &output);
        return array_size(nodes) - 1;
    }

    // Reorganize the primitive's nodes in order to be able to split them
    int primitive_in_b1 = 0;

    int mid = count0[minCostIndex] + start;
    for(int i = start; i < mid; i ++) {
        double3 offset = bbox_p_offset(bounds, primitives_info[i].center);
        int b = BUCKET_COUNT * offset.s[dim];
        if(b == BUCKET_COUNT)
            b = BUCKET_COUNT - 1;

        if(b > minCostIndex) {
            PrimitiveInfo tmp = primitives_info[i];
            Triangle tmp_tri = tris[i];

            primitives_info[i] = primitives_info[primitive_in_b1 + mid];
            tris[i] = tris[primitive_in_b1 + mid];

            primitives_info[primitive_in_b1 + mid] = tmp;
            tris[primitive_in_b1 + mid] = tmp_tri;

            primitive_in_b1 ++;
            i --;
        }
    }
    
    output = (BVHNode) {
        .bounds = bounds,
        .sons.x = buildBVHTreeSAH(obj, obj_meta, nodes, primitives_info, tris, start, mid),
        .sons.y = buildBVHTreeSAH(obj, obj_meta, nodes, primitives_info, tris, mid, end)
    };

    array_push(nodes, &output);

    int node_index = array_size(nodes) - 1;

    ((BVHNode*)nodes->array)[output.sons.x].parent = node_index;
    ((BVHNode*)nodes->array)[output.sons.y].parent = node_index;
    
    return node_index;
}

int buildBVHTreeFromObject(Object obj, ObjectMetadata *obj_meta)
{
    PrimitiveInfo node;
    Array nodes, primitives_info;

    nodes = array_create(sizeof(BVHNode));
    primitives_info = array_create(sizeof(PrimitiveInfo));

    for(uint i = 0; i < obj_meta->triangles_count; i++) {
        double3 points[3] = {
            obj.vertices[obj.triangles[i].vertices[0]],
            obj.vertices[obj.triangles[i].vertices[1]],
            obj.vertices[obj.triangles[i].vertices[2]]
        };

        node.bounds = (BBox3) {.min = points[0], .max = points[0]};
        node.bounds = bbox_p_union(node.bounds, points[1]);
        node.bounds = bbox_p_union(node.bounds, points[2]);
        node.center = bbox_center(node.bounds);

        array_push(&primitives_info, &node);
    }

    int root_id = buildBVHTreeSAH(obj, obj_meta, &nodes, primitives_info.array, obj.triangles, 0, obj_meta->triangles_count);

    obj_meta->tree = (BVHTree) {
        .nodes = nodes.array,
        .nodes_count = array_size(&nodes),
        .root = root_id
    };
    obj_meta->tree.nodes[root_id].parent = -1;

    array_free(&primitives_info);

    return 0;
}

uint BVHTree_depth(BVHNode* tree, int root)
{
    if(root == -1)
        return 0;
    
    return 1 + fmax(BVHTree_depth(tree, tree[root].sons.x), BVHTree_depth(tree, tree[root].sons.y));
}

void printObjectInfo(ObjectMetadata metadata)
{
    printf("================\n");
    printf("%ld triangles\n", metadata.triangles_count);
    printf("%ld vertices\n", metadata.vertices_count);
    printf("%ld normals\n", metadata.vertices_normal_count);
    printf("%ld texture coordonates\n", metadata.vertices_tex_count);
    printf("%ld meterials\n", metadata.materials_count);
    printf("BVH Tree depth: %ld\n", BVHTree_depth(metadata.tree.nodes, metadata.tree.root));
    printf("Bounds:\n");
    printf("- min: %lf %lf %lf\n",
        metadata.tree.nodes[metadata.tree.root].bounds.min.x,
        metadata.tree.nodes[metadata.tree.root].bounds.min.y,
        metadata.tree.nodes[metadata.tree.root].bounds.min.z);
    printf("- max: %lf %lf %lf\n",
        metadata.tree.nodes[metadata.tree.root].bounds.max.x,
        metadata.tree.nodes[metadata.tree.root].bounds.max.y,
        metadata.tree.nodes[metadata.tree.root].bounds.max.z);
    printf("Groups:\n");
    for(size_t i = 0; i < metadata.groups_count; i ++)
        printf(" - \"%s\"\n", metadata.groups_name[i]);
    printf("================\n");
}