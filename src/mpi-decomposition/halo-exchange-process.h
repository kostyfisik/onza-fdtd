#ifndef SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
#define SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
///
/// @file   halo-exchange-process.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Thu Apr 26 19:27:05 2012
///
/// @brief  Exchange borders of computational domain
#include <blitz/array.h>
#include "../simulation-core/basic-fdtd.h"
namespace onza {
  /// @brief Container for halo data to exchange.
  class HaloToExchange {
   public:
    /// @brief Start non-blocking communications.
    int StartNonBlockingExchange();
    /// @brief Finish exchange initiated with StartNonBlockingExchange().
    int FinishNonBlockingExchange();
   private:
    /// @brief My neighbours ranks
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions.
    /// @see HaloExchangeProcess::neighbours_ranks_
    int neighbours_ranks_[6];
    /// @brief Width of halo to exchange.
    int halo_width_;
  };  // end of class HaloToExchange
  /// @brief Class for MPI process.
  /// Contains computational domain borders data, methods to exchange
  /// and update borders.
  class HaloExchangeProcess {
   public:
    int Init();
    /// @breif Prepare simulation_core_ to start simulation.
    int InitSimulation();
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
    /// @brief Number of subdomains in each direction after domain
    /// decomposition.
    int64_t subdomains_[kDimensions];
    /// @brief Initialize HaloExchangeProcess::cartesian_grid_communicator_
    int StartCartesianGridCommunicator();
    /// @brief This process coords in Cartesian topology of MPI processes.
    int my_coords_[kDimensions];
    /// @brief My neighbours ranks
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions.
    int neighbours_ranks_[6];
    /// @brief Evaluate current process coords and neighbours ranks.
    ///
    /// Call after StartCartesianGridCommunicator()
    int EvaluateCoordsAndRanks();
    /// @brief Container for data to be exchange and methods to deal
    /// with this data.
    HaloToExchange halo_to_exchange_;
    // @}
    /// @name Computational section
    // @{
    /// @brief To carry out all simulation work.
    BasicSimulationCore simulation_core_;
    /// @brief Current process subdomain size in each dimension.
    ///
    /// size == number of vertices
    int64_t subdomain_size_[kDimensions];
    /// @breif Index of global model pointing to the beginnig of the
    /// subdomain.
    ///
    /// 0 index of subdomain should refer to same part of model
    /// as subdomain_start_index_
    int64_t subdomain_start_index_[kDimensions];
    /// @breif Index of global model pointing to the end of the
    /// subdomain.
    /// 
    /// @breif subdomain_size_-1 index of subdomain should refer to
    /// same part of model as subdomain_finish_index_
    int64_t subdomain_finish_index_[kDimensions];
    /// @brief Evaluate current process subdomain size.
    int EvaluateSubdomainSize();
    /// @brief Check subdomains start and finish indexes to be
    /// correlated with neighbours indexes.
    int CheckSubdomainIndexes();
    // @}
  };  // end of class HaloExchangeProcess
}  // end of namespace onza
#endif  // SRC_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
