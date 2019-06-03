/**
 * Graphics functions for controlling the animation of the demo.
 */
#pragma once

#include <SDL2/SDL.h>
#include "state.h"

class Graphics
{
public:
    /**
     * Default constructor. Uses default values for window size.
     */
    Graphics(void);

    /**
     * Parameterized constructor. Uses given width and height values for window size.
     */
    Graphics(int width, int height);

    /**
     * Initialize the graphics library. Return nonzero if it fails.
     */
    int init(void);

    /**
     * Cleanly exit the graphics library. Return nonzero if it fails.
     */
    int exit(void);

    /**
     * Updates the screen with the newest information. Returns nonzero if there is an error.
     */
    int update(State &state);

    /**
     * Returns whether we are done or not.
     */
    bool isdone(void);

private:
    bool done;
    int screen_width;
    int screen_height;

    SDL_Window *window;
    SDL_Surface *screensurface;
    SDL_Surface *background;

    /** Process the event queue, handling any user input. */
    int process_event_queue(void);

    /** Update the screen with drawings. */
    int update_image(State &state);
};
