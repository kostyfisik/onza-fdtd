///
/// @file   halo-exchange-process.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Thu Apr 26 19:13:58 2012
///
/// @brief  Exchange borders of computational domain
#include <mpi.h>
#include <cstdio>
#include "./halo-exchange-process.h"
#include "../common.h"
#include "../simulation-core/basic-fdtd.h"
namespace onza {
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::Init() {
    process_rank_ = MPI::COMM_WORLD.Get_rank();
    processes_total_number_ = MPI::COMM_WORLD.Get_size();
    printf("My MPI rank is %i from %i total\n",
           process_rank_, processes_total_number_);
    simulation_core_.simulation_input_config_.ReadConfig();
    int64_t length_x = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisX);
    int64_t length_y = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisY);
    int64_t length_z = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisZ);
    printf("HaloExchange length_x = %li, length_x = %li, length_x = %li\n",
           length_x, length_y, length_z);
    return 0;
  }  // end of HaloExchangeProcess::Init()
  /// @brief Run simulation domain decomposition
  ///
  /// With no information about cluster topology simply minimize surface
  /// of exchange border.
  /// @return 0 or error code.
  int  HaloExchangeProcess::RunDecomposition() {
    return 0;
  }  // end of HaloExchangeProcess::RunDecomposition
}  // end of namespace onza

