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
    simulation_input_config_.ReadConfig();
    int64_t length_x = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisX);
    int64_t length_y = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisY);
    int64_t length_z = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisZ);
    test_array_.resize(length_x, length_y, length_z);
    test_array_(0, 0, Range::all()) = 0;
    return 0;
  }  // end of BasicSimulationCore::Init
  /// @brief Read top level config file
  ///
  /// @todo Currently values to read from config file are hard coded
  /// in ReadConfig(). Read them from real config file.
  int SimulationInputConfig::ReadConfig() {
    // Length of whole model
    int64_t length_x = 10, length_y = 20, length_z = 30;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    status_ = static_cast<int>(kInputConfigAllDone);
    printf("Current status_ : %i \n", status_);  
    return 0;
  };  // end of SimulationInputConfig::ReadConfig
}  // end of namespace onza
