// Using Doxygen 1.8.0
///
/// @file   onza-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Wed Apr 25 15:01:02 2012
/// 
/// @brief  Top level calls for FDTD simulation, Doxygen mainpage description.
/// 
/// @mainpage Onza FDTD Documentation
///
/// [FDTD]: http://en.wikipedia.org/wiki/Finite-difference_time-domain_method "Go to Wikipedia"
/// [OpenMPI]: http://www.open-mpi.org
/// [Metamaterials Laboratory]: http://phoi.ifmo.ru/metamaterials
/// [NRU ITMO]: http://en.ifmo.ru "National Research University ITMO"
/// [Blitz++]: http://www.oonumerics.org/blitz
/// [CMake]: http://www.cmake.org "CMake"
/// [Ioffe Institute]: http://www.ioffe.ru/index_en.html "Ioffe Physical Technical Instute"
///
/// #Generic electromagnetic simulation software
///
/// Onza FDTD is a high performance electromagnetic simulation
/// software using [finite-difference time-domain (FDTD)
/// method][FDTD]. FDTD method generates very high computational
/// load. To get simulation results in reasonable time your should use
/// a powerful computer or a supercomputer cluster. Onza FDTD was designed
/// to run efficiently with:
///
///- multi-core processors
///- multiprocessor systems
///- clusters and supercomputers
///
/// Onza FDTD software is developed by [Metamaterials Laboratory] of
/// Photonics and Optical Informatics Department of [NRU ITMO]. It is
/// also supported by [Ioffe Institute].
///
/// Contact Ladutenko Konstantin <kostyfisik at gmail (.)  com>
/// with any questions about Onza FDTD.
///
/// #Prerequisites
///
///- MPI (to use Onza FDTD, to compile it)
///
///  Onza FDTD needs some MPI realization to be compiled and to be
///  run. It was developed using [OpenMPI]. For Debian/Ubuntu systems
///  it can be installed with:
///
///      # apt-get install openmpi-bin openmpi-doc libopenmpi-dev    
///
///- [Blitz++] library (to compile Onza FDTD)
///
///  For Debian/Ubuntu systems it can be installed with
///
///      # apt-get install libblitz0-dev
///
///- [CMake] build system (to compile Onza FDTD)
///
///  Use it on Linux/Unix and other operation systems.
#include <cstdio>
/// @todo blitz is internal implementation of model storge, move 
/// `#include <blitz/array.h>` there.
#include <blitz/array.h>
#include "mpi-decomposition/halo-exchange-process.h"
/// @brief Init, run and complete simulation.
///
/// @param argc 
/// @param argv well known input parameters
///
/// @return Zero by default.
int main(int argc, char *argv[])
{
  printf("Starting Onza!!!\n");
  HaloExchangeProcess halo_exchange_process;
  halo_exchange_process.init();
  return 0;
}
