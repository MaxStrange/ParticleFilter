#include <iostream>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "graphics.h"

/** Screen dimensions constants */
const int DEFAULT_SCREEN_WIDTH  = 640;
const int DEFAULT_SCREEN_HEIGHT = 480;

/** Colors */
const uint8_t ROBOT_RED     = 0x00;
const uint8_t ROBOT_GREEN   = 0x00;
const uint8_t ROBOT_BLUE    = 0xFF;
const uint8_t ROBOT_ALPHA   = 0xFF;

Graphics::Graphics(void)
{
    this->done = false;
    this->screen_width = DEFAULT_SCREEN_WIDTH;
    this->screen_height = DEFAULT_SCREEN_HEIGHT;
    this->window = nullptr;
    this->screensurface = nullptr;
    this->background = nullptr;
    this->renderer = nullptr;
}

Graphics::Graphics(int width, int height)
{
    this->done = false;
    this->screen_width = width;
    this->screen_height = height;
    this->window = nullptr;
    this->screensurface = nullptr;
    this->background = nullptr;
    this->renderer = nullptr;
}

int Graphics::get_screenheight(void) const
{
    return this->screen_height;
}

int Graphics::get_screenwidth(void) const
{
    return this->screen_width;
}

int Graphics::init(void)
{
    int err;

    err = SDL_Init(SDL_INIT_VIDEO);
    if (err)
    {
        err = __LINE__;
        goto fail;
    }

    this->window = SDL_CreateWindow("Particle Filter", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, this->screen_width, this->screen_height, SDL_WINDOW_SHOWN);
    if (this->window == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    this->renderer = SDL_CreateRenderer(this->window, -1, SDL_RENDERER_ACCELERATED);
    if (this->renderer == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    err = SDL_SetRenderDrawColor(this->renderer, 0xFF, 0xFF, 0xFF, 0xFF);
    if (err)
    {
        err = __LINE__;
        goto fail;
    }

    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG))
    {
        err = __LINE__;
        goto fail;
    }

fail:
    return err;
}

int Graphics::exit(void)
{
    SDL_DestroyRenderer(this->renderer);
    SDL_DestroyWindow(this->window);
    this->renderer = nullptr;
    this->window = nullptr;

    IMG_Quit();
    SDL_Quit();

    return 0;
}

int Graphics::update(const Robot &robot)
{
    int ret = 0;

    // Update based on state
    ret |= this->update_image(robot);

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

int Graphics::update_image(const Robot &robot)
{
    int err;
    SDL_Rect robotrect = { robot.get_xloc(), robot.get_yloc(), robot.get_width(), robot.get_height() };

    /* Clear the screen */
    err = SDL_SetRenderDrawColor(this->renderer, 0xFF, 0xFF, 0xFF, 0xFF);
    if (err)
    {
        goto fail;
    }

    err = SDL_RenderClear(this->renderer);
    if (err)
    {
        goto fail;
    }

    /* Paint the robot */
    SDL_SetRenderDrawColor(this->renderer, ROBOT_RED, ROBOT_GREEN, ROBOT_BLUE, ROBOT_ALPHA);
    SDL_RenderFillRect(this->renderer, &robotrect);

    /* Paint the particles */
    // TODO

    /* Update the screen */
    SDL_RenderPresent(this->renderer);

fail:
    return err;
}
