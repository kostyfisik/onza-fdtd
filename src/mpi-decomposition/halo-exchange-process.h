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
    /// @brief Run decomposition of simulation domain.
    int RunDecomposition();
   private:
    /// @brief Optized decomposition for star topology of supercomputer.
    int StarTopologyDecomposition(const int64_t all_processes,
                                  const int64_t length_x,
                                  const int64_t length_y,
                                  const int64_t length_z,
                                  int64_t &best_nx,
                                  int64_t &best_ny,
                                  int64_t &best_nz);
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
