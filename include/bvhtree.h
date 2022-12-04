#ifndef _BVH_TREE_H_
#define _BVH_TREE_H_

#include "bbox.h"
#include "typedef.h"

typedef struct BVHNode BVHNode;
struct BVHNode
{
    BBox3 bounds;
    int2 sons;
    uint triangle_start;
    uint triangle_end;
};

typedef struct BVHTree BVHTree;
struct BVHTree
{
    BVHNode *nodes;
    uint nodes_count;
    uint root;
};

#endif