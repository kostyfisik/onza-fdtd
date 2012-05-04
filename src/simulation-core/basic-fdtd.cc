///
/// @file   basic-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed May  2 13:42:32 2012
///
/// @brief  Simulation data and stepping algorithm
#include <mpi.h>
#include <blitz/array.h>
#include <cstdio>
#include "./basic-fdtd.h"
#include "../common.h"
namespace onza {
  using blitz::Range;
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  ///
  int  BasicSimulationCore::Init() {
    test_array_(10, 10, 10);
    test_array_(0, 0, Range::all()) = 0;
    simulation_input_config_.ReadSimulationInputConfig();
    printf("BasicSimulationCore::init() \n");
    return 0;
  }  // end of BasicSimulationCore::Init
  /// @brief Read top level config file
  ///
  /// @todo Currently values to read from config file are hard coded
  /// in ReadSimulationInputConfig(). Read them from real config file.
  int SimulationInputConfig::ReadSimulationInputConfig() {
    // Length of whole model
    int64_t length_x = 10, length_y = 10, length_z = 10;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    return 0;
  };  // end of SimulationInputConfig::ReadSimulationInputConfig
}  // end of namespace onza
