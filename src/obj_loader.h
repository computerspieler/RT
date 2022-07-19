#ifndef _OBJ_LOADER_H_
#define _OBJ_LOADER_H_

#include <stdio.h>

#include "array.h"
#include "object.h"
#include "tree.h"

int loadFromMtl(FILE *f, Node *materials_tree, Array *material_name, Array *materials);
int loadFromObj(FILE *f, Object *obj, ObjectMetadata *metadata);
//int mergeObjects(Object *obj1, ObjectMetadata *obj1_meta, Object *obj2, ObjectMetadata *obj2_meta, Object *output, ObjectMetadata *output_meta);

void printObjectInfo(ObjectMetadata metadata);

int buildBVHTreeFromObject(Object obj, ObjectMetadata *obj_meta);

#endif