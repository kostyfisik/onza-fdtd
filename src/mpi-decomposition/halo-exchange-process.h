#ifndef SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
#define SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
///
/// @file   halo-exchange-process.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Thu Apr 26 19:27:05 2012
///
/// @brief  Exchange borders of computational domain
#include "../simulation-core/basic-fdtd.h"
namespace onza {
  /// @brief Class for MPI process.
  /// Contains computational domain borders data, methods to exchange
  /// and update borders.
  class HaloExchangeProcess {
   public:
    int Init();
   private:
    /// @brief To carry out all non-communication work.
    BasicSimulationCore simulation_core_;
    /// @brief Fast access (without MPI call) to process rank of
    /// containing object. Should be set in init()
    int process_rank_;
    /// @brief Fast access (without MPI call) to total number of
    /// processes executing simulation. Should be set in init().
    int processes_total_number_;
  };  // end of class HaloExchangeProcess
}  // end of namespace onza
#endif  // SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
