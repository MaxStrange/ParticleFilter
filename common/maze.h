#pragma once

/** An x or y unit in the maze. Could be meters, feet, or whatever. */
typedef double maze_unit_t;

/**
 * The maze that the robot must navigate, using the particle filter.
 * Is initialized by a bitmap in the executable's directory.
 */
class Maze
{
public:
    Maze(void);
    ~Maze(void);

private:
};
