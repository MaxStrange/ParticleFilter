/**
 * Graphics functions for controlling the animation of the demo.
 */
#pragma once

#ifdef _WIN32
    #ifndef DEBUG
        #include "../SDL/include/SDL.h"
    #endif
#else
    #ifndef DEBUG
        #include <SDL2/SDL.h>
    #endif
#endif
#include "particlefilter.h"
#include "robot.h"

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

    int get_screenwidth(void) const;
    int get_screenheight(void) const;

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
    int update(const Robot &robot, const ParticleFilter &pf);

    /**
     * Returns whether we are done or not.
     */
    bool isdone(void);

private:
    bool done;
    int screen_width;
    int screen_height;

#ifndef DEBUG
    SDL_Window *window;
    SDL_Surface *screensurface;
    SDL_Surface *background;
    SDL_Renderer *renderer;
#endif

    /** Process the event queue, handling any user input. */
    int process_event_queue(void);

    /** Update the screen with drawings. */
    int update_image(const Robot &robot, const ParticleFilter &pf);
};
