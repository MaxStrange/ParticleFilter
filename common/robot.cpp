#include <iostream>
#include <random>
#include "robot.h"

#define MIN_VELOCITY_CHANGE -5
#define MAX_VELOCITY_CHANGE 5

Robot::Robot(void)
{
    this->xloc = 0;
    this->yloc = 0;
    this->xvel = 0;
    this->yvel = 0;
    this->height = 10;
    this->width = 10;

    std::random_device dev;
    this->rng = std::mt19937(dev());
    this->distribution = std::uniform_int_distribution<std::mt19937::result_type>(MIN_VELOCITY_CHANGE, MAX_VELOCITY_CHANGE);
}

Robot::~Robot(void)
{
}

int Robot::get_height(void) const
{
    return this->height;
}

int Robot::get_width(void) const
{
    return this->width;
}

int Robot::get_xloc(void) const
{
    return this->xloc;
}

int Robot::get_yloc(void) const
{
    return this->yloc;
}

void Robot::get_xy_estimate(int *x, int *y)
{
    // Return noisy estimate of my location
    *x = (this->distribution(this->rng) * 2) + this->get_xloc();
    *y = (this->distribution(this->rng) * 2) + this->get_yloc();
}

void Robot::get_v_estimate(int *vx, int *vy)
{
    // Return the noisy estimate of my velocity
    *vx = this->distribution(this->rng);
    *vy = this->distribution(this->rng);
}

void Robot::update(int height_limit, int width_limit)
{
    this->xloc += this->xvel;
    this->yloc += this->yvel;

    this->xvel = this->distribution(this->rng);
    this->yvel = this->distribution(this->rng);

    this->xloc = (this->xloc > width_limit) ? width_limit : this->xloc;
    this->xloc = (this->xloc < 0) ? 0 : this->xloc;
    this->yloc = (this->yloc > height_limit) ? height_limit : this->yloc;
    this->yloc = (this->yloc < 0) ? 0 : this->yloc;
}
