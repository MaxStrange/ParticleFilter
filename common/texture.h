#pragma once

#include <string>
#include <SDL2/SDL.h>

/**
 * This class wraps an SDL texture, providing lots of convenience around it.
 * Stolen from the lazyfoo tutorial.
 */
class Texture
{
public:
    Texture(void);
    ~Texture(void);

    void set_color(uint8_t r, uint8_t g, uint8_t b);
    void set_blendmode(SDL_BlendMode blending);
    void set_alpha(uint8_t alpha);

    /** Attempts to load the texture from the bitmap at the given location. */
    int load_from_file(const std::string &fpath, SDL_Renderer *renderer);

    /** Renders the texture at the given point in the screen. */
    int render(SDL_Renderer *renderer, int x, int y, const SDL_Rect *clip=nullptr, double angle=0.0, const SDL_Point *center=nullptr, SDL_RendererFlip flip=SDL_FLIP_NONE);

private:
    SDL_Texture *texture;

    int width;
    int height;

    /** Deallocates the texture */
    void free(void);
};
