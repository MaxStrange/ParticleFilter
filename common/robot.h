#pragma once

#include <random>

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

    int get_radius(void) const;
    int get_xloc(void) const;
    int get_yloc(void) const;

    /**
     * Update the robot's position by xvel and yvel and update
     * xvel and yvel by a small random amount. Assures the robot
     * does not go out of the limits provided.
     */
    void update(int height_limit, int width_limit);

private:
    // location
    int xloc;
    int yloc;

    // velocity
    int xvel;
    int yvel;

    // display
    int radius;

    std::uniform_int_distribution<std::mt19937::result_type> distribution;
    std::mt19937 rng;
};
