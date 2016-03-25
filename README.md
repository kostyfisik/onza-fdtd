Onza FDTD is a high performance electromagnetic simulation software
designed to run efficiently in parallel on multicore processors,
multiprocessor systems, clusters and supercomputers.

It only solves some predefined cases of Maxwell equation (usually anisotropic case) for 1,2,3D case an scales on HPC cluster.



Onza FDTD needs MPI installed to compile and run. 

Before use install Blitz++ library, MPI and Cmake. 

For Debian/Ubuntu systems single line install with

    apt-get install libblitz0-dev openmpi-bin openmpi-doc libopenmpi-dev cmake

./go.sh normaly should compile Onza and run FDTD simulation, edit it
to run Onza on your number of processes.
