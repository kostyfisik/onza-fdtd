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
    /// @breif accesor 
    int process_rank() {return process_rank_;}
   private:
    /// @brief Optized decomposition for star topology of supercomputer
    /// for 3D simulation model.
    int SimpleTopologyDecomposition3D(const int64_t all_processes,
                                  const int64_t in_length[],
                                  int64_t subdomains[]);
    /// @brief To carry out all non-communication work.
    BasicSimulationCore simulation_core_;
    /// @name MPI section
    // @{
    /// @brief Fast access (without MPI call) to process rank of
    /// containing object. Should be set in init()
    int process_rank_;
    /// @brief Fast access (without MPI call) to total number of
    /// processes executing simulation. Should be set in init().
    int processes_total_number_;
    /// @brief Cartesian grid communicator. Initiallized during model
    /// decomposition with RunDecomposition().
    MPI_Comm cartesian_grid_communicator_;
    /// @brief Initialize HaloExchangeProcess::cartesian_grid_communicator_ 
    int StartCartesianGridCommunicator(int64_t subdomains[]);
    /// @brief This process coords in Cartesian topology of MPI processes.
    int my_coords_[kDimensions];
    /// @brief My neighbors ranks
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions.
    int neighbors_ranks_[6];
    /// @brief Evaluate current process coords and neighbors ranks.
    ///
    /// Call after StartCartesianGridCommunicator()
    int EvaluateCoordsAndRanks();
    // @}
  };  // end of class HaloExchangeProcess
}  // end of namespace onza
#endif  // SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
