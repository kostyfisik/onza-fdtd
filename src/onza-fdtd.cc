// Using Doxygen 1.8.0 (with Markdown)
///
/// @file   onza-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed Apr 25 15:01:02 2012
///
/// @brief  Top level calls for FDTD simulation, Doxygen mainpage description.
#include <mpi.h>
#include <cstdio>
#include "mpi-decomposition/halo-exchange-process.h"
/// @brief Init, run and complete simulation.
///
/// @param argc
/// @param argv well known input parameters.
///
/// @return Zero by default.
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  onza::HaloExchangeProcess halo_exchange_process;
  int done_status = onza::kDone;
  while (1) {  // use break to report error with done_status.
    done_status = halo_exchange_process.Init();
    if (done_status != onza::kDone) break;
    done_status = halo_exchange_process.RunDecomposition();
    if (done_status != onza::kDone) break;
    done_status = halo_exchange_process.InitSimulation();
    break;
  }  // end of while breaked with errors
  MPI_Finalize();
  return done_status;
}  // end of main
// ************************************************************************* //
/// @page TopLevelAlgorithm Top level algorithm steps
///- Start exchange process (onza::HaloExchangeProcess::Init()). For each
/// MPI process:
/// * Get MPI runtime parameters
/// * Read simulation global input config
///  (onza::SimulationInputConfig::ReadConfig()).
/// * Cary out domain decomposition
///  onza::HaloExchangeProcess::RunDecomposition().  For early stages
///  of development we will assume star topology, decomposition is
///  optimized do minize volume of exchange data.
///- Start simulation on subdomain with help of onza::BasicSimulationCore
/// * Note: division could be several times longer than
///  multiplication. This way do not use division while time
///  stepping. All divisions should be done at inititalizing step.
/// @todo2 (For all project) Should be some logging system in Onza.
/// @todo3 (For all project) Perform out of range check for all ceil and float
/// operations.
// ************************************************************************* //
/// @page ChangeLog
/// #ChangeLog
/// @Section Version 0.0.1
///- CMake files are configured for ease of use.
// ************************************************************************* //
/// @page "Coding Style"
/// [Google style]: http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml
/// #Coding Style
/// @section Main
/// Use [Google style]!
///
/// @section Additional
///- Allways end class (methods, functions, namespaces, etc.) declaration
///  with comment
///
///        class GridInputConfig {
///          ...
///        };  // end of class GridInputConfig
///        int  BasicSimulationCore::Init() {
///          ...
///        }  // end of BasicSimulationCore::Init
///- All classes, containing input parameters for simulation should end
/// with *InputConfig*
// ************************************************************************* //
// /// \page page2 New page
// /// \tableofcontents
// /// Leading text.
// /// \section sec An example section
// /// This page contains the subsections \ref subsection1 and \ref subsection2.
// /// For more info see page \ref page2.
// /// \subsection subsection1 The first subsection
// /// Text.
// /// \subsection subsection2 The second subsection
// /// More text.
// ///
// ************************************************************************* //
/// @namespace onza
/// @brief Top level namespace of Onza FDTD project.
///
/// Also using names from MPI, blitz and std namespaces.
// ************************************************************************* //
// ************************************************************************* //
// ************************************************************************* //
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
/// #Parallel generic electromagnetic simulation software
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
///   Onza FDTD needs some MPI realization to be compiled and to be
///   run. It was developed using [OpenMPI]. For Debian/Ubuntu systems
///   it can be installed with:
///
///       # apt-get install openmpi-bin openmpi-doc libopenmpi-dev
///
///- [Blitz++] library (to compile Onza FDTD)
///
///   For Debian/Ubuntu systems it can be installed with
///
///       # apt-get install libblitz0-dev
///
///- [CMake] build system (to compile Onza FDTD)
///
///   Use it on Linux/Unix and other operation systems.
///
