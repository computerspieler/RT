#include <SDL/SDL_events.h>
#include <SDL/SDL_keysym.h>
#include <SDL/SDL_video.h>
#include <bits/types/struct_timeval.h>
#include <sys/time.h>

#include "cl_handler.h"
#include "camera.h"
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

const cl_double speed = 0.125f;

long long get_current_time()
{
	struct timeval te;
	gettimeofday(&te, NULL);
	return te.tv_sec * 1000LL + te.tv_usec / 1000;
}

int updateWindow(OpenCL_GeneralContext *cl_gen, OpenCL_ProgramContext *cl_prg, Camera *c)
{
	cl_int err;
	SDL_Event ev;
	int running = 1;

    size_t origin[3] = {0};
    size_t region[3] = {SCREEN_WIDTH, SCREEN_HEIGHT, 1};

	long long time_spent = 0;
	long long time_spent_opencl = 0;
	int nb_cycles = 0;
	long long date_last_refresh = get_current_time();

    TRY(opencl_add_output_image, exit,
        cl_gen, cl_prg, (cl_uchar4**) &screen->pixels, SCREEN_WIDTH, SCREEN_HEIGHT
    );
    
	while(running)
	{
		SDL_LockSurface(screen);
		long long start = get_current_time();
		opencl_run(cl_gen, cl_prg, origin, region);
		long long program_end = get_current_time();
		SDL_UnlockSurface(screen);
		
		//blitBuffer(screen, output, 0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
		SDL_Flip(screen);
		long long end = get_current_time();

		nb_cycles ++;
		time_spent += (end - start);
		time_spent_opencl += (program_end - start);
		
		if(end - date_last_refresh >= 1000L) {
			date_last_refresh = end;
			if(nb_cycles * 1000LL / time_spent > 0)
				printf("FPS: %lld\n", nb_cycles * 1000LL / time_spent);
			else
				printf("SPF: %lld\n", time_spent / 1000LL);
			printf("OpenCL: %lld%%\n", time_spent_opencl * 100LL / time_spent);
			time_spent = 0;
			time_spent_opencl = 0;
			nb_cycles = 0;
		}

		if(!SDL_PollEvent(&ev))
			continue;

		switch(ev.type)
		{
			case SDL_QUIT:
				running = 0;
				break;
			
			case SDL_KEYDOWN:
				switch(ev.key.keysym.sym) {
				case SDLK_z:			c->pos.z += speed; break;
				case SDLK_s:			c->pos.z -= speed; break;
				case SDLK_q:			c->pos.x -= speed; break;
				case SDLK_d:			c->pos.x += speed; break;
				case SDLK_SPACE:		c->pos.y -= speed; break;
				case SDLK_BACKSPACE:	c->pos.y += speed; break;
				
				default:;
				}
				break;
		}
	}

exit:
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

void blitBuffer(SDL_Surface* s, cl_uchar4* buffer, int x, int y, int w, int h)
{
	SDL_LockSurface(s);
	uint32_t *pixels = (uint32_t*) s->pixels;
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w; i ++){
            int k = i + j * w;
            pixels[x + i + (y + j) * s->w] =
				SDL_MapRGB(s->format, buffer[k].x, buffer[k].y, buffer[k].z);
        }

    SDL_UnlockSurface(s);
}