#pragma once

#include "maze.h"
#include "robot.h"
#include "particlefilter.h"

/**
 * This class contains all the state that the demonstration needs:
 * a maze, a robot, and a particle swarm. This is really meant
 * to just be a struct of things that are not really related, for
 * ease of passing all the things around and memory manipulation.
 */
struct State
{
public:
    /** This class is not responsible for the allocation/deallocation of pf */
    State(const ParticleFilter *const pf);
    ~State(void);

    Maze *maze;
    Robot *robot;
    const ParticleFilter *const pf;
};
