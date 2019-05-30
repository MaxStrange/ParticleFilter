/**
 * Much of this has been stolen most graciously from lazyfoo's SDL tutorial:
 * https://lazyfoo.net/tutorials/SDL/
 */
#include <SDL2/SDL.h>
#include <iostream>

/** Screen dimensions constants */
const int SCREEN_WIDTH  = 640;
const int SCREEN_HEIGHT = 480;

/** Exit macro: If not predicate, print errormsg and exit. */
#define CHECK(predicate, errormsg) do \
    {\
        if (!(predicate))\
        {\
            std::cout << errormsg << "Error: " << SDL_GetError() << "\n";\
            SDL_Quit();\
            exit(-1);\
        }\
    } while(0)

int main(void)
{
    auto err = SDL_Init(SDL_INIT_VIDEO);
    CHECK(err >= 0, "SDL could not initialize!");

    SDL_Window *window = SDL_CreateWindow("Particle Filter (Unoptimized)", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    CHECK(window != nullptr, "Window could not be created!");

    SDL_Surface *screensurface = SDL_GetWindowSurface(window);
    CHECK(screensurface != nullptr, "Could not get screen surface.");

    // Fill with white
    SDL_FillRect(screensurface, nullptr, SDL_MapRGB(screensurface->format, 0xFF, 0xFF, 0xFF));

    // Update the surface
    SDL_UpdateWindowSurface(window);

    // Wait two seconds
    SDL_Delay(2000);

    SDL_DestroyWindow(window);

    SDL_Quit();

    return 0;
}
