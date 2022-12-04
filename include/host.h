#ifndef _HOST_H_
#define _HOST_H_

#include "scene.h"
#include "camera.h"

typedef struct HostContext HostContext;
struct HostContext
{
    Camera camera;
    SceneMetadata scene_meta;
    int seed;
    int sampleCount;
};


#endif
