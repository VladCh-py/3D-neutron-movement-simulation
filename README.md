The 3D Monte-carlo Neutrons movement simulation.py simulates neutron transport in a 3D domain composed of different media.

Each neutron is initialized with a random starting position and direction within the source region. The source region is defined at runtime.
For each neutron, the free path length is sampled from an exponential distribution based on the macroscopic cross-sections of the medium. The neutron travels in a straight line until an interaction occurs, after which one of two events takes place: absorption by the medium (the neutron disappears) or scattering at a random angle relative to its previous direction.

During program execution:

*if enabled, neutron motion can be visualized in real time;
*simultaneously, in a separate window, a histogram is plotted representing the neutron flux as a function of distance r from the cylinder center.
The plots are generated using the Matplotlib library.

Neutron flux calculation:
The neutron flux is computed as the number of neutrons crossing a unit area per unit time (or, in this case, within a unit radial interval in cylindrical coordinates).
During each simulation cycle, the current positions of all neutrons are recorded and used to build a histogram of their spatial distribution.
The height of each histogram bin is proportional to the neutron density in that region of space, which effectively represents the local neutron flux as a function of the radial distance from the axis (in this particular case, since cylindrical geometry is used).

Geometry and input data:
In this specific example, a cylindrical geometry is considered, defined via an external OBJ file.
This file contains the description of the region boundaries (media and source). The OBJ file can be created in any modeling software, but it must not contain non-closed surfaces, otherwise the simulation may become unstable.

After each simulation step, absorbed (“dead”) neutrons are replaced with new ones to maintain a constant total number of particles in the system, ensuring the stationarity of the process.



Folder "pictures" contains screenshots of the program, the 3D models used, and some of the results.

source_line.obj
The files "source_line.obj" and "materials.mtl" contain the geometry and description of the materials of the cylindrical medium. These files are used in the example (must be placed in the same folder with the program).
