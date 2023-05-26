#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "bbox.h"
#include "bvhtree.h"
#include "vec3.h"
#include "typedef.h"
#include "material.h"
#include "scene_loader.h"
#include "array.h"
#include "scene.h"
#include "tree.h"
#include "light.h"

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
            
            mat = (Material) {0};
            fscanf(f, "%s %d\n", buffer, (int*) &mat.type);
            str_length = strlen(buffer);
            name = calloc(str_length + 1, sizeof(char));
            memcpy(name, buffer, str_length * sizeof(char));
            name[str_length] = 0;
            array_push(material_name, &name);
            add_word(materials_tree, name, array_size(materials));

            if(!initialized)
                initialized = true;
        }
        else {
            switch(mat.type) {
                case MATERIAL_NONE: break;
                case MATERIAL_LAMBERTIAN:
                    if(strcmp(buffer, "R") == 0)
                        fscanf(f, FLOAT_FMT "\n", &mat.l.R);
                    break;
                
                case MATERIAL_GLASS:
                    if(strcmp(buffer, "IOR") == 0)
                        fscanf(f, FLOAT_FMT "\n", &mat.g.IOR);
                    break;

                default:                    
                    if(strcmp(buffer, "Sr") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.specular_roughness);
                    }
                    else if(strcmp(buffer, "SIOR") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.specular_IOR);
                    }
                    else if(strcmp(buffer, "M") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.metalness);
                    }
                    else if(strcmp(buffer, "M") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.metalness);
                    }
                    else if(strcmp(buffer, "T") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.transmission);
                    }
                    else if(strcmp(buffer, "Td") == 0) {
                        fscanf(f, FLOAT_FMT "\n", &mat.c.transmission_dispersion);
                    }
                    break;
            }
        }
    }
    
    if(initialized)
        array_push(materials, &mat);

    return 0;
}

void subdivideTriangles(Scene *obj, SceneMetadata *metadata, int subdivide)
{
	if(subdivide <= 0)	return;

	obj->vertices = realloc(obj->vertices, (metadata->vertices_count + metadata->triangles_count) * sizeof(vec3));

	for(size_t i = 0; i < metadata->triangles_count; i ++) {
	    vec3 p0 = obj->vertices[obj->triangles[i].vertices[0]];
	    vec3 p1 = obj->vertices[obj->triangles[i].vertices[1]];
	    vec3 p2 = obj->vertices[obj->triangles[i].vertices[2]];

		obj->vertices[i + metadata->vertices_count].x = (p0.x + p1.x + p2.x) / 3.;
		obj->vertices[i + metadata->vertices_count].y = (p0.y + p1.y + p2.y) / 3.;
		obj->vertices[i + metadata->vertices_count].z = (p0.z + p1.z + p2.z) / 3.;
	}

	obj->vertices_tex = realloc(obj->vertices_tex, (metadata->vertices_tex_count + metadata->triangles_count) * sizeof(vec2));

	for(size_t i = 0; i < metadata->triangles_count; i ++) {
	    vec2 p0 = obj->vertices_tex[obj->triangles[i].uv[0]];
	    vec2 p1 = obj->vertices_tex[obj->triangles[i].uv[1]];
	    vec2 p2 = obj->vertices_tex[obj->triangles[i].uv[2]];

		obj->vertices_tex[i + metadata->vertices_tex_count].x = (p0.x + p1.x + p2.x) / 3.;
		obj->vertices_tex[i + metadata->vertices_tex_count].y = (p0.y + p1.y + p2.y) / 3.;
	}


	obj->triangles = realloc(obj->triangles, 3 * metadata->triangles_count * sizeof(Triangle));

	for(size_t i = 0; i < metadata->triangles_count; i ++) {
		obj->triangles[metadata->triangles_count + 2 * i] = obj->triangles[i];
		obj->triangles[metadata->triangles_count + 2 * i + 1] = obj->triangles[i];

		obj->triangles[i].uv[2] = i + metadata->vertices_tex_count;
		obj->triangles[i].vertices[2] = i + metadata->vertices_count;

		obj->triangles[metadata->triangles_count + 2 * i].uv[0] = i + metadata->vertices_tex_count;
		obj->triangles[metadata->triangles_count + 2 * i].vertices[0] = i + metadata->vertices_count;
		
		obj->triangles[metadata->triangles_count + 2 * i + 1].uv[1] = i + metadata->vertices_tex_count;
		obj->triangles[metadata->triangles_count + 2 * i + 1].vertices[1] = i + metadata->vertices_count;
	}


	metadata->vertices_count += metadata->triangles_count;
	metadata->vertices_tex_count += metadata->triangles_count;
	metadata->triangles_count *= 3;

	subdivideTriangles(obj, metadata, subdivide - 1);
}

int loadFromObj(FILE *f, Scene *obj, SceneMetadata *metadata, bool use_uv)
{
    FILE *fm;
    vec3 vertex;
	Light light;
    vec2 vertex_uv;
    Triangle triangle;

    int res;
    char buffer[1024];
    char *name;
    uint current_material = -1;
    uint str_length = 0;
    Array vertices, tex, triangles, groups, material_name, materials, lights;
    Node *materials_tree, *n;

    assert(obj);
    assert(metadata);

    vertices = array_create(sizeof(vec3));
    tex = array_create(sizeof(vec2));
    groups = array_create(sizeof(char*));
    triangles = array_create(sizeof(Triangle));
    materials = array_create(sizeof(Material));
    materials_tree = malloc(sizeof(Node));
    bzero(materials_tree, sizeof(Node));
    material_name = array_create(sizeof(char*));
    material_name = array_create(sizeof(char*));
    lights = array_create(sizeof(Light));

    if(!f) {
        perror("File");
        return -1;
    }

    while(true) {
        res = fscanf(f, "%s", buffer);
        if(res == EOF)
            break;

        if(strcmp(buffer, "v") == 0) {
            fscanf(f, FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n", &vertex.x, &vertex.y, &vertex.z);
            array_push(&vertices, &vertex);
        } else if(strcmp(buffer, "vt") == 0) {
            fscanf(f, FLOAT_FMT " " FLOAT_FMT "\n", &vertex_uv.x, &vertex_uv.y);
            array_push(&tex, &vertex_uv);
        } else if(strcmp(buffer, "f") == 0) {
            res = fscanf(f, "%ld/%ld %ld/%ld %ld/%ld\n",
                &triangle.vertices[0], &triangle.uv[0],
                &triangle.vertices[1], &triangle.uv[1],
                &triangle.vertices[2], &triangle.uv[2]);

            triangle.material = current_material;
            triangle.group = array_size(&groups) - 1;
            
            triangle.vertices[0] --;
            triangle.vertices[1] --;
            triangle.vertices[2] --;
            
            triangle.uv[0] --;
            triangle.uv[1] --;
            triangle.uv[2] --;

			triangle.original_id = array_size(&triangles);

            if(res != 6) {
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

            res = loadFromMtl(fm, materials_tree, &material_name, &materials);
            fclose(fm);

            if(res)
                goto err;
        } else if(strcmp(buffer, "lp") == 0) {
            fscanf(f, FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "\n", &light.pos.x, &light.pos.y, &light.pos.z, &light.I);
            light.type = LIGHT_POINT;
            array_push(&lights, &light);
        } else if(strcmp(buffer, "la") == 0) {
            fscanf(f,
				FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT " "
				FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT " "
				FLOAT_FMT " " FLOAT_FMT " "
				FLOAT_FMT "\n",
				&light.pos.x, &light.pos.y, &light.pos.z,
				&light.area_n.x, &light.area_n.y, &light.area_n.z,
				&light.area_width, &light.area_height,
				&light.I);
            light.type = LIGHT_AREA;

			vec3_normalize(light.area_n);
			vec3_build_coordonate_system(light.area_n, &light.area_tan, &light.area_bitan);
            array_push(&lights, &light);
        } else
            while (fgetc(f) != '\n');
    }

    free_tree(materials_tree);

    *metadata = (SceneMetadata) {
        .groups_count = array_size(&groups),
        .groups_name = groups.array,
        .materials_count = array_size(&materials),
        .triangles_count = array_size(&triangles),
        .vertices_tex_count = array_size(&tex),
        .vertices_count = array_size(&vertices),
        .use_uv = use_uv,
		.lights_count = array_size(&lights)
    };

    *obj = (Scene) {
        .triangles = triangles.array,
        .vertices = vertices.array,
        .vertices_tex = tex.array,
        .materials = materials.array,
		.lights = lights.array
    };

	subdivideTriangles(obj, metadata, SUBDIVIDE_COUNT);

    buildBVHTreeFromObject(*obj, metadata);

    for(size_t i = 0; i < metadata->triangles_count; i ++) {
	    vec3 p0 = obj->vertices[obj->triangles[i].vertices[0]];
	    vec3 p1 = obj->vertices[obj->triangles[i].vertices[1]];
	    vec3 p2 = obj->vertices[obj->triangles[i].vertices[2]];

		assert(obj->triangles[i].uv[0] >= 0);
		assert(obj->triangles[i].uv[1] >= 0);
		assert(obj->triangles[i].uv[2] >= 0);
		assert(obj->triangles[i].uv[0] < (long int) metadata->vertices_tex_count);
		assert(obj->triangles[i].uv[1] < (long int) metadata->vertices_tex_count);
		assert(obj->triangles[i].uv[2] < (long int) metadata->vertices_tex_count);

	    vec2 uv0 = obj->vertices_tex[obj->triangles[i].uv[0]];
	    vec2 uv1 = obj->vertices_tex[obj->triangles[i].uv[1]];
	    vec2 uv2 = obj->vertices_tex[obj->triangles[i].uv[2]];

        vec3 dp02 = vec3_diff(p0, p2);
        vec3 dp12 = vec3_diff(p1, p2);

        vec2 duv02 = VEC2(
            uv0.x - uv2.x,
            uv0.y - uv2.y
        );
        vec2 duv12 = VEC2(
            uv1.x - uv2.x,
            uv1.y - uv2.y
        );

	    obj->triangles[i].normal = vec3_cross(vec3_diff(p2, p0), vec3_diff(p1, p0));

		obj->triangles[i].area = vec3_norm(obj->triangles[i].normal) * .5;
		obj->triangles[i].normal = vec3_normalize(obj->triangles[i].normal);

        Float det = duv02.x * duv12.y - duv02.y * duv12.x;

        if(det == 0)
            vec3_build_coordonate_system(obj->triangles[i].normal, &obj->triangles[i].dpdu, &obj->triangles[i].dpdv);
        else {
            obj->triangles[i].dpdu = vec3_normalize(
                vec3_diff(vec3_smul(duv12.y, dp02), vec3_smul(duv02.y, dp12))
            );
            obj->triangles[i].dpdv = vec3_normalize(
                vec3_diff(vec3_smul(duv02.x, dp12), vec3_smul(duv12.x, dp02))
            );
        }
    }

	printf("==== BVH Tree ====\n");
	printf("Node count: %d\n", metadata->tree.nodes_count);
	printf("Tree count: %ld\n", metadata->tree.nodes_count * sizeof(BVHNode));

    printf("Scene file sucessfully loaded\n");

    return 0;

err:
    free_tree(materials_tree);
    array_free(&vertices);
    array_free(&tex);
    array_free(&triangles);
    array_free(&groups);
    array_free(&materials);
    array_free(&lights);

    return -1;
}

typedef struct PrimitiveInfo PrimitiveInfo;
struct PrimitiveInfo
{
    BBox3 bounds;
    vec3 center;
};

int buildBVHTreeSAH(Scene obj, SceneMetadata *scene_meta, Array *nodes, PrimitiveInfo *primitives_info, Triangle *tris, int start, int end)
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
        vec3 offset = bbox_p_offset(bounds, primitives_info[i].center);
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
    Float cost[BUCKET_COUNT - 1];
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
            cost[i] = 0.125 + (Float) (count0[i] * bbox_surface_area(b0) + count1 * bbox_surface_area(b1)) / bbox_surface_area(bounds);
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
        vec3 offset = bbox_p_offset(bounds, primitives_info[i].center);
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
        .sons.x = buildBVHTreeSAH(obj, scene_meta, nodes, primitives_info, tris, start, mid),
        .sons.y = buildBVHTreeSAH(obj, scene_meta, nodes, primitives_info, tris, mid, end)
    };

    array_push(nodes, &output);

    int node_index = array_size(nodes) - 1;
    
    return node_index;
}

int buildBVHTreeFromObject(Scene obj, SceneMetadata *scene_meta)
{
    PrimitiveInfo node;
    Array nodes, primitives_info;

    nodes = array_create(sizeof(BVHNode));
    primitives_info = array_create(sizeof(PrimitiveInfo));

    for(uint i = 0; i < scene_meta->triangles_count; i++) {
        vec3 points[3] = {
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

    int root_id = buildBVHTreeSAH(obj, scene_meta, &nodes, primitives_info.array, obj.triangles, 0, scene_meta->triangles_count);

    scene_meta->tree = (BVHTree) {
        .nodes = nodes.array,
        .nodes_count = array_size(&nodes),
        .root = root_id
    };

    array_free(&primitives_info);

    return 0;
}

uint BVHTree_depth(BVHNode* tree, int root)
{
    if(root == -1)
        return 0;
    
    return 1 + fmax(BVHTree_depth(tree, tree[root].sons.x), BVHTree_depth(tree, tree[root].sons.y));
}

void printBVHTree(BVHNode* tree, int node, int depth)
{
    if(node == -1) return;

    for(int i = 0; i < depth; i ++)
        putchar(' ');

    printf("ID: %d; Bounds: [Min: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "; ",
        node,
        tree[node].bounds.min.x,
        tree[node].bounds.min.y,
        tree[node].bounds.min.z);
    printf("Max: " FLOAT_FMT " " FLOAT_FMT " " FLOAT_FMT "]; ",
        tree[node].bounds.max.x,
        tree[node].bounds.max.y,
        tree[node].bounds.max.z);
    
    printf("Triangles: [Start = %d; End = %d]\n",
        tree[node].triangle_start, tree[node].triangle_end);

    printBVHTree(tree, tree[node].sons.x, depth + 1);
    printBVHTree(tree, tree[node].sons.y, depth + 1);
}

void printObjectInfo(SceneMetadata metadata)
{
    printf("================\n");
    printf("%ld triangles\n", metadata.triangles_count);
    printf("%ld sommets\n", metadata.vertices_count);
    printf("%ld coordonnées de textures\n", metadata.vertices_tex_count);
    printf("%ld matériaux\n", metadata.materials_count);
    printf("Profondeur de l'arbre BVH: %u\n", BVHTree_depth(metadata.tree.nodes, metadata.tree.root));
    //printBVHTree(metadata.tree.nodes, metadata.tree.root, 0);
    printf("Groupes:\n");
    for(size_t i = 0; i < metadata.groups_count; i ++)
        printf(" - \"%s\"\n", metadata.groups_name[i]);
    printf("================\n");
}

void scene_delete(Scene *obj, SceneMetadata *metadata)
{
	assert(obj);
	assert(metadata);

	free(obj->vertices);
	free(obj->vertices_tex);
	free(obj->triangles);
	free(obj->materials);
	free(obj->lights);
	free(metadata->groups_name);
	free(metadata->tree.nodes);
}