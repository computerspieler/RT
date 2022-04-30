#include "render.h"

SDL_Surface* screen;

int initWindow()
{
	if(SDL_Init(SDL_INIT_VIDEO) < 0)
		SDL_FATAL(SDL_Init);

	screen = SDL_SetVideoMode(
			SCREEN_WIDTH,
			SCREEN_HEIGHT,
			32, SDL_HWSURFACE
		);

	if(!screen)
		SDL_FATAL(SDL_SetVideoMode);

	SDL_WM_SetCaption("RayTracing", NULL);

	return 0;
}

int updateWindow(char* output)
{
	SDL_Event ev;
	int running = 1;

    blitBuffer(screen, output, 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

    SDL_Flip(screen);

	while(running)
	{
		if(!SDL_PollEvent(&ev))
			continue;

		switch(ev.type)
		{
			case SDL_QUIT:
				running = 0;
				break;
		}
	}

	return 0;
}

void fillRect(SDL_Surface* s, int x, int y, int w, int h, Color c)
{
	for(int j = y; j < y+h; j++)
		for(int i = x; i < x+w; i++)
			drawPixel(s, i, j, c);		
}

void drawPixel(SDL_Surface *s, int x, int y, Color c)
{
	uint32_t *pixels = (uint32_t*) s->pixels;
	pixels[x + y * s->w] = SDL_MapRGB(s->format, c.r, c.g, c.b);
}

void blitBuffer(SDL_Surface* s, char* buffer, int x, int y, int w, int h)
{
	SDL_LockSurface(s);
	uint32_t *pixels = (uint32_t*) s->pixels;
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w; i ++){
            int k = i + j * w;
            pixels[x + i + (y + j) * s->w] =
				SDL_MapRGB(s->format, buffer[4*k+1], buffer[4*k+2], buffer[4*k+3]);
        }

    SDL_UnlockSurface(s);
}