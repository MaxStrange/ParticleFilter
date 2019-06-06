#include <string>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "texture.h"

Texture::Texture(void)
{
    this->height = 0;
    this-> width = 0;
    this->texture = nullptr;
}

Texture::~Texture(void)
{
    this->free();
}

void Texture::set_color(uint8_t r, uint8_t g, uint8_t b)
{
    SDL_SetTextureColorMod(this->texture, r, g, b);
}

void Texture::set_blendmode(SDL_BlendMode blending)
{
    SDL_SetTextureBlendMode(this->texture, blending);
}

void Texture::set_alpha(uint8_t alpha)
{
    SDL_SetTextureAlphaMod(this->texture, alpha);
}

int Texture::load_from_file(const std::string &fpath, SDL_Renderer *renderer)
{
    int err;
    SDL_Texture *newtexture = nullptr;

    // Override this Texture's current contents
    this->free();

    SDL_Surface *loadedsurface = IMG_Load(fpath.c_str());
    if (loadedsurface == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    err = SDL_SetColorKey(loadedsurface, SDL_TRUE, SDL_MapRGB(loadedsurface->format, 0, 0xFF, 0xFF));
    if (err)
    {
        goto fail;
    }

    newtexture = SDL_CreateTextureFromSurface(renderer, loadedsurface);
    if (newtexture == nullptr)
    {
        err = __LINE__;
        goto fail;
    }

    this->width = loadedsurface->w;
    this->height = loadedsurface->h;

    SDL_FreeSurface(loadedsurface);

fail:
    return err;
}

int Texture::render(SDL_Renderer *renderer, int x, int y, const SDL_Rect *clip/*=nullptr*/, double angle/*=0.0*/, const SDL_Point *center/*=nullptr*/, SDL_RendererFlip flip/*=SDL_FLIP_NONE*/)
{
    SDL_Rect renderquad = {
        x,
        y,
        this->width,
        this->height,
    };

    if (clip != nullptr)
    {
        renderquad.w = clip->w;
        renderquad.h = clip->h;
    }

    return SDL_RenderCopyEx(renderer, this->texture, clip, &renderquad, angle, center, flip);
}

///////////////////////////////////////////////////////////
//////////////    PRIVATE FUNCTIONS   /////////////////////
///////////////////////////////////////////////////////////

void Texture::free(void)
{
    if (this->texture != nullptr)
    {
        SDL_DestroyTexture(this->texture);
        this->texture = nullptr;
        this->width = 0;
        this->height = 0;
    }
}
