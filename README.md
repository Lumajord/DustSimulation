# Dust Simulation
CMake project; Updated version of the dust simulation code by Alex Seizinger. The code simulates dust agglomerates by computing the surface forces of dust monomers based on the model by [Wada et al. 2007](https://ui.adsabs.harvard.edu/abs/2007ApJ...661..320W/abstract).  
The code runs on CPU and GPU via Cuda, the CPU version has a GUI created with QT.

Two example simulations with the code: compression of a dust agglomerate inside a box and a catastrophic collision between two dust agglomerates. 
<img src="https://github.com/Lumajord/dustsimulation/blob/main/box.jpg" width="250"> <img src="https://github.com/Lumajord/dustsimulation/blob/main/Collision.png" width="350">
