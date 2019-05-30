#include <iostream>
#include <SDL2/SDL.h>
#include "graphics.h"

/** Screen dimensions constants */
const int DEFAULT_SCREEN_WIDTH  = 640;
const int DEFAULT_SCREEN_HEIGHT = 480;

Graphics::Graphics(void)
{
    this->done = false;
    this->screen_width = DEFAULT_SCREEN_WIDTH;
    this->screen_height = DEFAULT_SCREEN_HEIGHT;
}

Graphics::Graphics(int width, int height)
{
    this->done = false;
    this->screen_width = width;
    this->screen_height = height;
}

int Graphics::init(void)
{
    int err;
    SDL_Window *window = nullptr;
    SDL_Surface *screensurface = nullptr;

    err = SDL_Init(SDL_INIT_VIDEO);
    if (err)
    {
        goto fail;
    }

    window = SDL_CreateWindow("Particle Filter (Unoptimized)", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, this->screen_width, this->screen_height, SDL_WINDOW_SHOWN);
    if (window == nullptr)
    {
        goto fail;
    }

    screensurface = SDL_GetWindowSurface(window);
    if (screensurface == nullptr)
    {
        goto fail;
    }

fail:
    return err;
}

int Graphics::exit(void)
{
    SDL_Quit();

    return 0;
}

int Graphics::update(void)
{
    SDL_Event e;

    // Check if we have an event in the queue
    auto got_event = SDL_PollEvent(&e);

    // Process the event if we got one
    if (got_event)
    {
        switch (e.type)
        {
            case SDL_QUIT:
                this->done = true;
                break;
            default:
                break;
        }
    }

    return 0;
}

bool Graphics::isdone(void)
{
    return this->done;
}
