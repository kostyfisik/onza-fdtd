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
    simulation_core_.Init();
    return 0;
  }  // end of HaloExchangeProcess::Init()
}  // end of namespace onza

