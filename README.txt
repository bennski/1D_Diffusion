"README.txt"

File: 1D_Diffusion.py
Author: Ben K. D. Pearce
Created: July 4, 2016
Updated: November 11, 2018

***Prerequisite packages***
1. python 2.7
2. ffmpeg

Description:
-The following code was written to simulate the outflow of molecules from a meteorite sitting in a pond.
-The code produces an animation of the outflow.

Boundary conditions:
-Neumann at the center of the meteorite (r=0), no molecules can leave from this boundary.
-Open at the edge of the meteorite (r=rmax), molecules can leave from this boundary.

Initial condition:
-Constant concentration from r=0 to r=rmax, except for a smoothly decreasing concentration to zero at the outer edge of the meteorite (rmax).


