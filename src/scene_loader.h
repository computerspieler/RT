#ifndef _OBJ_LOADER_H_
#define _OBJ_LOADER_H_

#include <stdio.h>

#include "array.h"
#include "scene.h"
#include "bvhtree.h"
#include "tree.h"

int loadFromMtl(FILE *f, Node *materials_tree, Array *material_name, Array *materials);
int loadFromObj(FILE *f, Scene *obj, SceneMetadata *metadata, bool use_uv);

void scene_delete(Scene *obj, SceneMetadata *metadata);

void printObjectInfo(SceneMetadata metadata);
uint BVHTree_depth(BVHNode* tree, int root);

int buildBVHTreeFromObject(Scene obj, SceneMetadata *scene_meta);

#endif