#ifndef _HOST_H_
#define _HOST_H_

#include "scene.h"
#include "camera.h"

typedef struct HostContext HostContext;
struct HostContext
{
    Camera camera;
	int2 mouse_pos;
    SceneMetadata scene_meta;
    uint seed;
    int lambda;
	uint sampleCount;
	bool simpleView;
	Float max_threshold;
};


#endif
