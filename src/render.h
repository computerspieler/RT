#ifndef _RENDER_H_
#define _RENDER_H_

#include <CL/cl.h>
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <SDL/SDL.h>
#include <SDL/SDL_video.h>

#include "cl_handler.h"
#include "camera.h"

typedef struct
{
	int r;
	int g;
	int b;
} Color;

#define SDL_FATAL(func)						\
{											\
	printf(#func ": %s\n", SDL_GetError());	\
	exit(EXIT_FAILURE);						\
}

#define SCREEN_WIDTH	120
#define SCREEN_HEIGHT	120

int initWindow();
int updateWindow(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, Camera *c);
void fillRect(SDL_Surface*, int x, int y, int w, int h, Color);
void drawPixel(SDL_Surface*, int x, int y, Color c);
void blitBuffer(SDL_Surface*, cl_uchar4*, int x, int y, int w, int h);

#endif
