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
    this->radius = 10;

    std::random_device dev;
    this->rng = std::mt19937(dev());
    this->distribution = std::uniform_int_distribution<std::mt19937::result_type>(MIN_VELOCITY_CHANGE, MAX_VELOCITY_CHANGE);
}

Robot::~Robot(void)
{
}

int Robot::get_radius(void) const
{
    return this->radius;
}

int Robot::get_xloc(void) const
{
    return this->xloc;
}

int Robot::get_yloc(void) const
{
    return this->yloc;
}

void Robot::update(int height_limit, int width_limit)
{
    this->xloc += this->xvel;
    this->yloc += this->yvel;

    this->xvel = this->distribution(rng);
    this->yvel = this->distribution(rng);

    this->xloc = (this->xloc > width_limit) ? width_limit : this->xloc;
    this->xloc = (this->xloc < 0) ? 0 : this->xloc;
    this->yloc = (this->yloc > height_limit) ? width_limit : this->yloc;
    this->yloc = (this->yloc < 0) ? 0 : this->yloc;

    std::cout << "Position: (" << this->xloc << ", " << this->yloc << ")" << std::endl;
}
