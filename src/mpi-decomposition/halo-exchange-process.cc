///
/// @file   halo-exchange-process.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Thu Apr 26 19:13:58 2012
///
/// @brief  Exchange borders of computational domain
#include <mpi.h>
#include <cstdio>
#include <cmath>
#include "./halo-exchange-process.h"
#include "../common.h"
#include "../simulation-core/basic-fdtd.h"
namespace onza {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Evaluate current process coords and neighbors ranks.
  int HaloExchangeProcess::EvaluateCoordsAndRanks() {
    MPI_Cart_coords(cartesian_grid_communicator_, process_rank_,
                    kDimensions, my_coords_);
    int upwards_shift = 1, downwards_shift = -1;
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      int rank_source, rank_dest;
      MPI_Cart_shift(cartesian_grid_communicator_, axis,
                     upwards_shift, &rank_source, &rank_dest);
      neighbors_ranks_[axis] = rank_dest;
      MPI_Cart_shift(cartesian_grid_communicator_, axis,
                     downwards_shift, &rank_source, &rank_dest);
      neighbors_ranks_[axis+3] = rank_dest;
    }
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Run decomposition of simulation domain.
  ///
  /// With no information about cluster topology simply minimize surface
  /// of exchange border using StarTopologyDecomposition3D().
  ///
  ///  @todo1 Domain decomposition should be correlated with topology
  ///  of cluster, starting from intra node cores\processors layout
  ///  with nested levels of borders exchange.
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::RunDecomposition() {
    if (simulation_core_.simulation_input_config_.
        status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;
    int64_t length[kDimensions];
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      length[axis] = simulation_core_.simulation_input_config_.
        grid_input_config_.get_total_grid_length(static_cast<Axis>(axis));
    }  // end of for axis
    const int64_t all_processes = processes_total_number_;
    int64_t subdomains[kDimensions];
    SimpleTopologyDecomposition3D(all_processes, length, subdomains);    
    if (StartCartesianGridCommunicator(subdomains) != kDone)
      return kErrorProcessNotInGrid;
    if (process_rank_==0) {
      printf("Total MPI processes = %i, used = %li, ",
             processes_total_number_, MultiplyComponents(subdomains));
      printf("subdomains (%li,%li,%li), domain length (%li,%li,%li)\n",
             subdomains[kAxisX],
             subdomains[kAxisY], subdomains[kAxisZ],
             length[kAxisX],
             length[kAxisY], length[kAxisZ]);
    }
    EvaluateCoordsAndRanks();
    /// @todo3 Remove output in HaloExchangeProcess::RunDecomposition()
    // printf("coords (%2i,%2i,%2i) : process_rank_: %2i    X(%2i,%2i) Y(%2i,%2i) Z(%2i,%2i)\n",
    //        my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ],
    //        process_rank_,
    //        neighbors_ranks_[kBorderRight], neighbors_ranks_[kBorderLeft],
    //        neighbors_ranks_[kBorderTop], neighbors_ranks_[kBorderBottom],
    //        neighbors_ranks_[kBorderFront], neighbors_ranks_[kBorderBack]);
    return kDone;
  }  // end of HaloExchangeProcess::RunDecomposition
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize cartesian grid communicator
  ///
  /// @param[in] subdomains[] Count of subdomains in each direction.
  /// @see HaloExchangeProcess::StarTopologyDecomposition3D() 
  /// @return 0 or error code.
  int HaloExchangeProcess::StartCartesianGridCommunicator(int64_t subdomains[]) {
    int dimensions_number = kDimensions;
    int processes_in_each_dimension[kDimensions];
    int is_grid_periodic[kDimensions];
    for (int i = kAxisX; i < kDimensions; ++i) {
      processes_in_each_dimension[i] = subdomains[i];
      is_grid_periodic[i] = 0; // no periodicity
    }
    bool reorder_processes_ranks = true;
    // Set cartesian_grid_communicator_
    MPI_Cart_create(MPI_COMM_WORLD, dimensions_number,
                    processes_in_each_dimension,
                    is_grid_periodic, reorder_processes_ranks,
                    &cartesian_grid_communicator_);
    if (cartesian_grid_communicator_ == MPI_COMM_NULL)
      return kErrorProcessNotInGrid;
    MPI_Comm_rank(cartesian_grid_communicator_, &process_rank_);
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::Init() {
    MPI_Comm_size(MPI_COMM_WORLD, &processes_total_number_);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    int init_status = simulation_core_.Init();
    if (init_status != kDone) return init_status;
    return kDone;
  }  // end of HaloExchangeProcess::Init()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Optized decomposition for star topology of supercomputer
  /// for 3D simulation model.
  ///
  /// Try to make computational volume for each MPI process to be cubic.
  /// For small number of processes (<10) we are trying to utilize
  /// them all (even if single domain become far from to be cubic).
  /// Algorithm is consant time, it directly compares fixed number of guessed
  /// 'good' decompositions.
  ///
  /// @param[in] all_processes Number of MPI processes available. Number
  /// of subdomains should be as close as possible to it.
  /// @param[in] in_length[] Size of volume for decomposition in each dimension
  /// @param[out] subdomains[] Count of subdomains in each direction.
  /// Suppozed to be the best ones for selected optimization.
  int HaloExchangeProcess::SimpleTopologyDecomposition3D(
      const int64_t all_processes,
      const int64_t in_length[],
      int64_t subdomains[]) {
    // arrays redefinement to instantate the template of MultiplyComponents()
    int64_t length[kDimensions], best_n[kDimensions];
    double normalized_length[kDimensions];
    int64_t floor_normalized_length[kDimensions];
    int64_t ceil_normalized_length[kDimensions];
    for (int i = kAxisX; i < kDimensions; ++i) {
        length[i] = in_length[i];
        normalized_length[i] = 0;
        floor_normalized_length[i] = 0;
        ceil_normalized_length[i] = 0;
      }
    // Prepare for guessing best decomposition
    double single_cell_optimal_volume = static_cast<double>
      (MultiplyComponents(length)) / static_cast<double>(all_processes);
    int optimal_dimenstion = kDimensions;
    double optimal_length = pow(single_cell_optimal_volume, 1.0/optimal_dimenstion);
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (length[i]/optimal_length < 1) optimal_dimenstion -= 1;
      normalized_length[i] = 1;
    }
    optimal_length = pow(single_cell_optimal_volume, 1.0/optimal_dimenstion);
    // Set floor and ceil normalized_length guess.
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (floor_normalized_length[i] == 1 &&
          ceil_normalized_length[i] == 1) {
        continue;  // floor and ceil values were evaluated earlier.
      }
      normalized_length[i] = length[i]/optimal_length;
      floor_normalized_length[i]
        = static_cast<int64_t>(floor(normalized_length[i]));
      // should be at least one subdomain in every valid direction.
      if (floor_normalized_length[i] < 1) floor_normalized_length[i] = 1;
      // for volumes too long in one direction
      // floor value can be very big
      if (floor_normalized_length[i] > all_processes)
        floor_normalized_length[i] = all_processes;
      ceil_normalized_length[i]
        = static_cast<int64_t>(ceil(normalized_length[i]));
      // nuber of subdomains in every direction should not
      // exceed total number of MPI processes
      if (ceil_normalized_length[i] > all_processes)
        ceil_normalized_length[i] = all_processes;
    }  // end for noralized length calculation
    // Selecting "guessed" values for decomposition,
    // 2 integer values are less than the optimal float value,
    // 3 integer values are bigger than the optimal value.
    const int variant_size = 5;
    // number of processes (subdomains) in selected dimension.
    int64_t variant[kDimensions][variant_size];
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (floor_normalized_length[i] == 1 &&
          ceil_normalized_length[i] == 1) {
        for (int j = 0; j < variant_size; ++j) 
          variant[i][j] = 1;
        best_n[i] = 1;
        continue;
      }
      variant[i][0] = floor_normalized_length[i]-1 < 1 ?
        1 : floor_normalized_length[i]-1;
      variant[i][1] = floor_normalized_length[i];
      variant[i][2] = ceil_normalized_length[i];
      variant[i][3] = ceil_normalized_length[i]+1 > all_processes ?
        all_processes : ceil_normalized_length[i]+1;
      variant[i][4] = ceil_normalized_length[i]+2 > all_processes ?
        all_processes : ceil_normalized_length[i]+2;
      // special case - big number of processes and small length
      for (int j = 0; j < variant_size; ++j)
        variant[i][j] = variant[i][j] > length[i] ? length[i] : variant[i][j];
      best_n[i] = variant[i][0];
    }
    // Total number of subdomains for current step of optimization
    // Trying to make it as close as possible to all_processes value.
    int64_t best_size = MultiplyComponents(best_n);
    // Total surface between subdomains for current step. The less is better.
    int64_t best_norm
      = (best_n[kAxisX]-1) * length[kAxisY] * length[kAxisZ]
      + (best_n[kAxisY]-1) * length[kAxisX] * length[kAxisZ]
      + (best_n[kAxisZ]-1) * length[kAxisX] * length[kAxisY];
    // Selecting best decomposition from the guessed ones.
    for (int i = 0; i < variant_size; ++i) {
      for (int j = 0; j < variant_size; ++j) {
        for (int k = 0; k < variant_size; ++k) {
          // Calculate determinating values for new decomposition
          int64_t new_size
            = variant[kAxisX][i] * variant[kAxisY][j] * variant[kAxisZ][k];
          int64_t new_norm
            = (variant[kAxisX][i]-1) * length[kAxisY] * length[kAxisZ]
            + (variant[kAxisY][j]-1) * length[kAxisX] * length[kAxisZ]
            + (variant[kAxisZ][k]-1) * length[kAxisX] * length[kAxisY];
          if (new_size > all_processes) continue;
          // Select new decomposition to be best if it is really better
          if (new_size == best_size && new_norm < best_norm) {
            best_norm = new_norm;
            best_n[kAxisX] = variant[kAxisX][i];
            best_n[kAxisY] = variant[kAxisY][j];
            best_n[kAxisZ] = variant[kAxisZ][k];
          }  // end of if new best norm and same size
          if (new_size > best_size) {
            best_size = new_size;
            best_norm = new_norm;
            best_n[kAxisX] = variant[kAxisX][i];
            best_n[kAxisY] = variant[kAxisY][j];
            best_n[kAxisZ] = variant[kAxisZ][k];
          }  // end of if new best size
        }  // end of for k
      }  // end of for j
    }  // end of for i
    for (int i = kAxisX; i < kDimensions; ++i) subdomains[i] = best_n[i];
    return kDone;
  }  // end of HaloExchangeProcess::StarTopologyDecomposition
}  // end of namespace onza

