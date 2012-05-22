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
  /// @brief Evaluate current process subdomain size.
  ///
  /// Take into account information about number of subdomains in
  /// each dimensions, with of PML layer, it`s relative
  /// computational difficulty.
  int HaloExchangeProcess::EvaluateSubdomainSize() {
    // Set aliases.
    int pml_width = simulation_core_.simulation_input_config_.pml_width();
    double pml_computational_ratio = simulation_core_.simulation_input_config_.
      pml_computational_ratio();
    int boundary_condition[kDimensions * 2];
    for (int border = kBorderLeft; border < kDimensions*2; ++border) 
      boundary_condition[border] = simulation_core_.simulation_input_config_.
        boundary_condition(static_cast<BorderPosition>(border));
    // Split for each axis.
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      // Number of vertices.
      double total_size = simulation_core_.simulation_input_config_.
        grid_input_config_.get_total_grid_length(static_cast<Axis>(axis));
      // Number of equivalent vertices.
      double total_computational_size = total_size;
      // For axis with PML the computational size is bigger.
      // PML is part of total grid.
      if (boundary_condition[axis] == kBoudaryConditionPML)
        total_computational_size += pml_width*(pml_computational_ratio - 1.0);
      if (boundary_condition[axis+kDimensions] == kBoudaryConditionPML)
        total_computational_size += pml_width*(pml_computational_ratio - 1.0);
      // @todo3 check if subdomain_computational_size is not too small.
      // Optimal number of equivalent vertices per subdomain dimension.
      double subdomain_computational_size =
        total_computational_size/static_cast<double>(subdomains_[axis]);
      // Evaluate grid coords from current MPI process coords.
      double subdomain_start_coord = -1;
      double subdomain_finish_coord = -1;
      // For wide PML/too many processes exist pure PML subdomains.
      double pml_subdomain_size = subdomain_computational_size
        /pml_computational_ratio;
      int64_t  pml_subdomains_number_from_left = 0;
      int64_t  pml_subdomains_number_from_right = 0;
      double transient_subdomain_pml_vertices_from_left = 0;
      double transient_subdomain_pml_vertices_from_right = 0;
      // Starting from first (left, bottom, back) side of domain.
      if (boundary_condition[axis] == kBoudaryConditionPML) {
        pml_subdomains_number_from_left
          = floor(pml_width/pml_subdomain_size);
        transient_subdomain_pml_vertices_from_left = pml_width -
          pml_subdomain_size * pml_subdomains_number_from_left;
      } // end of if first boundary is PML
      // Second (right, top, front) side of domain.
      if (boundary_condition[axis+kDimensions] == kBoudaryConditionPML) {
        pml_subdomains_number_from_right
          = pml_subdomains_number_from_left;
        transient_subdomain_pml_vertices_from_right =
          transient_subdomain_pml_vertices_from_left;
      } // end of if second boundary is PML
      // Inside PML - pure pml subdomains.
      /// Overlaping PML in same directions should be forbidden
      /// in config file with SimulationInputConfig::ReadConfig().
      if (my_coords_[axis] < pml_subdomains_number_from_left) {
        subdomain_start_coord = pml_subdomain_size * my_coords_[axis];
        subdomain_finish_coord =
          pml_subdomain_size * (my_coords_[axis] + 1);
      }  // end of if pure pml near first boundary
      if (subdomains_[axis] - 1 - my_coords_[axis]
          < pml_subdomains_number_from_right) {
        subdomain_start_coord = total_size -
          pml_subdomain_size *
          (subdomains_[axis] - my_coords_[axis]);
        subdomain_finish_coord = total_size -
          pml_subdomain_size *
          (subdomains_[axis] - 1 - my_coords_[axis]);
      }  // end of if pure pml near second boundary
      // Transient subdomain, partly with PML.
      double left_transient_subdomain_start_coord =
        pml_subdomain_size * pml_subdomains_number_from_left;
      double right_transient_subdomain_finish_coord = total_size 
        - pml_subdomain_size * pml_subdomains_number_from_right;
      // Transient subdomains are bigger than pure pml subdomain and
      // less or equal to computational subdomain size.
      double left_transient_subdomain_finish_coord =
        left_transient_subdomain_start_coord + subdomain_computational_size - 
        transient_subdomain_pml_vertices_from_left*(pml_computational_ratio - 1.0);
      double right_transient_subdomain_start_coord =
        right_transient_subdomain_finish_coord - subdomain_computational_size + 
        transient_subdomain_pml_vertices_from_right*(pml_computational_ratio - 1.0);
      if (my_coords_[axis] == pml_subdomains_number_from_left) {
        subdomain_start_coord = left_transient_subdomain_start_coord;
        subdomain_finish_coord = left_transient_subdomain_finish_coord;
        // Special case: left and right PML are in same transient subdomain.
        if (subdomain_finish_coord > right_transient_subdomain_finish_coord)
          subdomain_finish_coord = right_transient_subdomain_finish_coord;
      }  // end of if transient subdomain near first boundary
      if (subdomains_[axis] - 1 - my_coords_[axis] == pml_subdomains_number_from_right
          // Exclude special case: left and right PML are in same transient subdomain.
          && subdomain_finish_coord < 0) {
        subdomain_finish_coord = right_transient_subdomain_finish_coord;
        subdomain_start_coord =  right_transient_subdomain_start_coord;
      }  // end of if transient subdomain near second boundary
      // Inner area.
      if (my_coords_[axis] > pml_subdomains_number_from_left &&
          subdomains_[axis] - 1 - my_coords_[axis]
          > pml_subdomains_number_from_right) {
        subdomain_start_coord = left_transient_subdomain_finish_coord
          + subdomain_computational_size *
            (my_coords_[axis] - pml_subdomains_number_from_left - 1);
        subdomain_finish_coord = subdomain_start_coord
          + subdomain_computational_size;
      }  // end of if subdomian in inner area (without PML)
      // Output.
      int zeros = 0;
      int non_zero_axis = -1;
      for (int i = 0; i < kDimensions; ++i)
        if (my_coords_[i] == 0) zeros += 1;
        else non_zero_axis = i;
      if ((zeros == 2 && axis == non_zero_axis) || zeros == 3) {
        printf("(%2i,%2i,%2i) axis[%i] computational: total vert %7.2f, sub vert%7.2f, \
start %7.2f finish %7.2f\n",
               my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ],
               axis, total_computational_size, subdomain_computational_size,
               ceil(subdomain_start_coord),
               ceil(subdomain_finish_coord) - 1);
      }  // end if process with coords output
      subdomain_size_[axis] = 0;
    }  // end of for axis
    return kDone;
  }  // end of HaloExchangeProcess::EvaluateSubdomainSize()
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
    // (0) Prepare.
    if (simulation_core_.simulation_input_config_.
        status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;
    int64_t length[kDimensions], orig_length[kDimensions];
    // Original length of simulation domain
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      orig_length[axis] = simulation_core_.simulation_input_config_.
        grid_input_config_.get_total_grid_length(static_cast<Axis>(axis));
      length[axis] = orig_length[axis];
    }
    int pml_width = simulation_core_.simulation_input_config_.pml_width();
    double pml_computational_ratio = simulation_core_.simulation_input_config_.
      pml_computational_ratio();
    int boundary_condition[kDimensions * 2];
    for (int border = kBorderLeft; border < kDimensions*2; ++border) 
      boundary_condition[border] = simulation_core_.simulation_input_config_.
        boundary_condition(static_cast<BorderPosition>(border));
    // Length of simulation domain with computational ajustment.
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      if (boundary_condition[axis] == kBoudaryConditionPML)
        length[axis] += ceil(pml_width*(pml_computational_ratio - 1.0));
      if (boundary_condition[axis+kDimensions] == kBoudaryConditionPML)
        length[axis] += ceil(pml_width*(pml_computational_ratio - 1.0));
    }  // end of for dimensions
    const int64_t all_processes = processes_total_number_;
    // (1) Topology decomposition.
    SimpleTopologyDecomposition3D(all_processes, length, subdomains_);
    if (process_rank_ == 0) {
      printf("Used MPI procs %li of %i (%.3g%%), ",
             MultiplyComponents(subdomains_), processes_total_number_,
             MultiplyComponents(subdomains_)/double(processes_total_number_)*100);
      printf("subdomains_ (%li,%li,%li), domain length (%li,%li,%li) equiv ->(%li,%li,%li)\n",
             subdomains_[kAxisX],
             subdomains_[kAxisY], subdomains_[kAxisZ],
             orig_length[kAxisX],
             orig_length[kAxisY], orig_length[kAxisZ],
             length[kAxisX],
             length[kAxisY], length[kAxisZ]);
    }  // end of if process_rank_ == 0 output
    // (2) Start new communicator.
    if (StartCartesianGridCommunicator() != kDone)
      return kErrorProcessNotInGrid;
    // (3) Find neighbors.
    EvaluateCoordsAndRanks();
    /// @todo3 Remove output in HaloExchangeProcess::RunDecomposition()
    // printf("(%2i,%2i,%2i) rank: %2i  X(%2i,%2i) Y(%2i,%2i) Z(%2i,%2i)\n",
    //        my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ],
    //        process_rank_,
    //        neighbors_ranks_[kBorderLeft], neighbors_ranks_[kBorderRight],
    //        neighbors_ranks_[kBorderBottom], neighbors_ranks_[kBorderTop],
    //        neighbors_ranks_[kBorderBack],  neighbors_ranks_[kBorderFront]);
    // (4) Find current process subdomain size.
    EvaluateSubdomainSize();
    return kDone;
  }  // end of HaloExchangeProcess::RunDecomposition
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
      neighbors_ranks_[axis+kDimensions] = rank_dest;
    }  // end of for dimensions
    // Set periodical boundary conditions.
    int boundary_condition[kDimensions * 2];
    for (int border = kBorderLeft; border < kDimensions*2; ++border) 
      boundary_condition[border] = simulation_core_.simulation_input_config_.
        boundary_condition(static_cast<BorderPosition>(border));
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (boundary_condition[border] != kBoudaryConditionPeriodical) continue;
      int axis = border % kDimensions;
      int direction_to_border = border < kDimensions ? 1 : -1;
      int max_shift = (subdomains_[axis] - 1);
      int rank_source, rank_dest;
      if (my_coords_[axis] == 0 && direction_to_border < 0) {        
        MPI_Cart_shift(cartesian_grid_communicator_, axis,
                       max_shift, &rank_source, &rank_dest);
        neighbors_ranks_[border] = rank_dest;
        continue;
      }
      if (my_coords_[axis] == max_shift && direction_to_border > 0) {
        MPI_Cart_shift(cartesian_grid_communicator_, axis,
                       -max_shift, &rank_source, &rank_dest);
        neighbors_ranks_[border] = rank_dest;
        continue;
      }
    }  // end of for dimensions
    return kDone;
  }  // end of HaloExchangeProcess::EvaluateCoordsAndRanks()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize cartesian grid communicator
  ///
  /// @see HaloExchangeProcess::StarTopologyDecomposition3D()
  /// @return 0 or error code.
  int HaloExchangeProcess::StartCartesianGridCommunicator() {
    int dimensions_number = kDimensions;
    int processes_in_each_dimension[kDimensions];
    int is_grid_periodic[kDimensions];
    for (int i = kAxisX; i < kDimensions; ++i) {
      processes_in_each_dimension[i] = subdomains_[i];
      // Set no periodicity. If periodicity is needed it should be set
      // up inside EvaluateCoordsAndRanks()
      is_grid_periodic[i] = 0;  
    }  // end of for dimensions
    bool reorder_processes_ranks = true;
    // Set new MPI communicator: cartesian_grid_communicator_.
    MPI_Cart_create(MPI_COMM_WORLD, dimensions_number,
                    processes_in_each_dimension,
                    is_grid_periodic, reorder_processes_ranks,
                    &cartesian_grid_communicator_);
    if (cartesian_grid_communicator_ == MPI_COMM_NULL)
      return kErrorProcessNotInGrid;
    MPI_Comm_rank(cartesian_grid_communicator_, &process_rank_);
    return kDone;
  }  // end of HaloExchangeProcess::StartCartesianGridCommunicator()
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
    int optimal_dimension = kDimensions;
    double optimal_length = pow(single_cell_optimal_volume,
                                1.0/optimal_dimension);
    int is_dimension_reduced[kDimensions];
    single_cell_optimal_volume = 1;
    // Find out dimenstions to reduce using full 3D decomposition
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (length[i]/optimal_length < 1) {
        optimal_dimension -= 1;
        is_dimension_reduced[i] = 1;
        continue;
      }
      is_dimension_reduced[i] = 0;
      single_cell_optimal_volume *= length[i];
    }
    // Re-evaluate optimization parameters for reduced (if present) case
    single_cell_optimal_volume /= static_cast<double>(all_processes);
    optimal_length = pow(single_cell_optimal_volume, 1.0/optimal_dimension);
    /// %todo3 Remove debug output
    // if (process_rank_ == 0) {
    //   printf("optimal_length %.3g \n", optimal_length);
    //   printf("single_cell_optimal_volume %.3g \n", single_cell_optimal_volume);
    //   printf("optimal_dimension %i \n", optimal_dimension);      
    // }
    // Set floor and ceil normalized_length guess.
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (is_dimension_reduced[i] == 1) {  // Reduced dimenstion is not optimized.
        floor_normalized_length[i] = 1;
        ceil_normalized_length[i] = 1;
        continue;  
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
    // 4 integer values are less than the optimal float value,
    // 3 integer values are bigger than the optimal value.
    const int variant_size = 7; 
    // Propable number of processes (subdomains) in selected dimension.
    int64_t variant[kDimensions][variant_size];
    for (int i = kAxisX; i < kDimensions; ++i) {
      if (is_dimension_reduced[i] == 1) {
        // For reduced dimensions there is no variants.
        for (int j = 0; j < variant_size; ++j)
          variant[i][j] = 1;
        best_n[i] = 1;
        continue;
      }  // end of if dimension is reduced
      variant[i][0] = floor_normalized_length[i]-3 < 1 ?
        1 : floor_normalized_length[i]-3;
      variant[i][1] = floor_normalized_length[i]-2 < 1 ?
        1 : floor_normalized_length[i]-2;
      variant[i][2] = floor_normalized_length[i]-1 < 1 ?
        1 : floor_normalized_length[i]-1;
      variant[i][3] = floor_normalized_length[i];
      variant[i][4] = ceil_normalized_length[i];
      variant[i][5] = ceil_normalized_length[i]+1 > all_processes ?
        all_processes : ceil_normalized_length[i]+1;
      variant[i][6] = ceil_normalized_length[i]+2 > all_processes ?
        all_processes : ceil_normalized_length[i]+2;
      // Special case - big number of processes and small length.
      for (int j = 0; j < variant_size; ++j)
        variant[i][j] = variant[i][j] > length[i] ? length[i] : variant[i][j];
      best_n[i] = variant[i][0];
    }  // end of for all dimensions
    // Total number of subdomains for current step of optimization.
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
          // Select new decomposition to be best if it is really better.
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
