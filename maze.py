"""
Create a "maze.bmp" based on data written here in human-readable format.

Maze stolen from particle_filter_demo (see README).
"""
import imageio
import numpy as np
import skimage

# 0 - empty square
# 1 - occupied square
# 2 - occupied square with a beacon at each corner, detectable by the robot

maze_data = ( ( 1, 1, 0, 0, 2, 0, 0, 0, 0, 1 ),
              ( 1, 2, 0, 0, 1, 1, 0, 0, 0, 0 ),
              ( 0, 1, 1, 0, 0, 0, 0, 1, 0, 1 ),
              ( 0, 0, 0, 0, 1, 0, 0, 1, 1, 2 ),
              ( 1, 1, 0, 1, 1, 2, 0, 0, 1, 0 ),
              ( 1, 1, 1, 0, 1, 1, 1, 0, 2, 0 ),
              ( 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 ),
              ( 1, 2, 0, 1, 1, 1, 1, 0, 0, 0 ),
              ( 0, 0, 0, 0, 1, 0, 0, 0, 1, 0 ),
              ( 0, 0, 1, 0, 0, 2, 1, 1, 1, 0 ))

if __name__ == "__main__":
    imagedata = np.array(maze_data) * 100.0
    imagedata = skimage.transform.resize(imagedata, (500, 500), anti_aliasing=False)
    imagedata[imagedata >= 125.0] = 200.0
    imagedata[(imagedata >= 50.0) & (imagedata < 125.0)] = 100.0
    imagedata[imagedata < 50.0] = 0.0
    imagedata = imagedata.astype(np.uint8)
    imageio.imwrite("maze.bmp", imagedata)
