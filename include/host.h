#ifndef _HOST_H_
#define _HOST_H_

#include "scene.h"

typedef struct HostContext HostContext;
struct HostContext
{
	int2 mouse_pos;
    SceneMetadata scene_meta;
    uint seed;
};


#endif
