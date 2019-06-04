#pragma once

#include "maze.h"

/**
 * The Robot is the agent that moves around the maze with the help of the particle filter.
 * It has sensors that it uses to get some noisy information about the environment,
 * and actuators, which it uses to move in a random direction.
 */
class Robot
{
public:
    Robot(void);
    ~Robot(void);

private:
    maze_unit_t xloc;
    maze_unit_t yloc;
    maze_unit_t xvel;
    maze_unit_t yvel;
};
