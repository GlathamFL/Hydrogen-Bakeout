# Hydrogen-Bakeout
Matlab code that is used to determine hydrogen concentration in various geometries made of 316 stainless steel.

Run in either Octave or MATLAB.

You can use this code to pick a geometry and assign dimensions for spheres, plates, or tubes. You can then choose a tempeature and time to perform a bakeout for.
The code will output a graph of the hydrogen concentration throughout a line in the solid, as well as the amount of moles of hydrogen left in the solid.
The idea is to minimize the amount of hydrogen in the solid so that the part does not degass hydrogen into a process like MOCVD or PVD.

The code can be adapated to other materials and diffusing atoms as long as you can the material's parameters related to Fick's laws.

Some assumptions made for the plate:
There is no diffusion of hydrogen out of the side pertaining to the thickness of the material.
This means hydrogen only diffuses out of (or into) the faces exposed to the process side and non-process side.
This code is also a numerical solution, and it not analytical in nature. This may mean some error is present in the output.

Some assumptions made for the tube:
The formula was taken from Crank's text on the mathematics of diffusion.
It assumes the initial concentration of hydrogen throughout the solid is constant.
It also assumes there is no diffusion out of the top and bottom of the solid.

Some assumptions made for the sphere:
The formula was taken from Crank's text on the mathematics of diffusion.
It assumes the initial concentration of hydrogen throughout the solid is constant.
