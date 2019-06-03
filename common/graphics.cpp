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
    this->window = nullptr;
    this->screensurface = nullptr;
    this->background = nullptr;
}

Graphics::Graphics(int width, int height)
{
    this->done = false;
    this->screen_width = width;
    this->screen_height = height;
    this->window = nullptr;
    this->screensurface = nullptr;
    this->background = nullptr;
}

int Graphics::init(void)
{
    int err;

    err = SDL_Init(SDL_INIT_VIDEO);
    if (err)
    {
        goto fail;
    }

    this->window = SDL_CreateWindow("Particle Filter (Unoptimized)", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, this->screen_width, this->screen_height, SDL_WINDOW_SHOWN);
    if (this->window == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    this->screensurface = SDL_GetWindowSurface(this->window);
    if (this->screensurface == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    this->background = SDL_LoadBMP("maze.bmp");
    if (this->background == nullptr)
    {
        err = __LINE__;
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

int Graphics::update(State &state)
{
    int ret = 0;

    // Update based on state
    ret |= this->update_image(state);

    // Update based on event queue
    ret |= this->process_event_queue();

    return ret;
}

bool Graphics::isdone(void)
{
    return this->done;
}

////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Private Functions ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

int Graphics::process_event_queue(void)
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

int Graphics::update_image(State &state)
{
    int err;

    // TODO: Update based on state

    err = SDL_BlitSurface(this->background, nullptr, this->screensurface, nullptr);
    if (err)
    {
        goto fail;
    }

    err = SDL_UpdateWindowSurface(this->window);
    if (err)
    {
        goto fail;
    }

fail:
    return err;
}
