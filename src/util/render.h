#ifndef _RENDER_H_
#define _RENDER_H_

#include <SDL/SDL.h>
#include <SDL/SDL_video.h>

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

void fillRect(SDL_Surface*, int x, int y, int w, int h, Color);
void drawPixel(SDL_Surface*, int x, int y, Color c);
void blitBuffer(SDL_Surface*, char*, int x, int y, int w, int h);

#endif
