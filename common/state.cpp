#include "state.h"

State::State(const ParticleFilter *const pf)
    : pf(pf)
{
    this->maze = new Maze();
    this->robot = new Robot();
}

State::~State(void)
{
    delete this->maze;
    delete this->robot;
}
