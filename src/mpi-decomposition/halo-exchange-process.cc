///
/// @file   halo-exchange-process.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of Onza FDTD.
///
/// Onza FDTD is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// Onza FDTD is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Onza FDTD.  If not, see <http://www.gnu.org/licenses/>.
/// @date   Thu Apr 26 19:13:58 2012
///
/// @brief  Exchange borders of computational domain
#include <mpi.h>
#include <blitz/array.h>
#include <cstdio>
#include <cmath>
#include <string>
#include "./halo-exchange-process.h"
#include "../common.h"
#include "../simulation-core/basic-fdtd.h"
#include "../profiling/timer.h"
namespace onza {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Start non-blocking communications.
  ///
  /// @warning  It is repeated MANY times!
  void HaloToExchange::StartNonBlockingExchange() {
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      int send_size = number_of_elements_to_send_[border];
      // receive border in current process was oppositly labeled in
      // sending process
      int direction_tag = border < kDimensions
        ? border + kDimensions : border - kDimensions;
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      // Case of border transfer from one side of subdomain to
      // opposite side is treated during preparing borders in
      // BasicSimulationCore::PrepareBordersToSend()
      if (neighbours_ranks_[border] == process_rank_) continue;
      MPI_Irecv(received_borders_[border], send_size, MPI_DOUBLE,
                neighbours_ranks_[border], direction_tag,
                cartesian_grid_communicator_, &irecv_request_[border]);
    }  // end of for borders receive
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      int send_size = number_of_elements_to_send_[border];
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      if (neighbours_ranks_[border] == process_rank_) continue;
      int direction_tag = border;
      MPI_Isend(borders_to_send_[border], send_size, MPI_DOUBLE,
                neighbours_ranks_[border], direction_tag,
                cartesian_grid_communicator_, &isend_request_[border]);
    }  // end of for borders send
  }  // end of HaloToExchange::StartNonBlockingExchange()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Finish exchange initiated with StartNonBlockingExchange().
  ///
  /// @warning  It is repeated MANY times!
  void HaloToExchange::FinishNonBlockingExchange() {
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      if (neighbours_ranks_[border] == process_rank_) continue;
      // For debug it could be helpfull to have some status parser...
      // MPI_Wait(&isend_request_[border], &status);
      // MPI_Wait(&irecv_request_[border], &status);
      MPI_Wait(&isend_request_[border], MPI_STATUS_IGNORE);
      MPI_Wait(&irecv_request_[border], MPI_STATUS_IGNORE);
    }  // end of for borders wait
  }  // end of HaloToExchange::FinishNonBlockingExchange()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Init exchange buffers parameters.
  ///
  /// Get addresses and sizes for input and output buffers to send and
  /// receive borders.
  int HaloToExchange::Init(
        blitz::Array<blitz::Array<double, 1+kDimensions>, 1> borders_to_send,
        blitz::Array<blitz::Array<double, 1+kDimensions>, 1> received_borders,
        int process_rank, int neighbours_ranks[],
        MPI_Comm &cartesian_grid_communicator) {
    process_rank_ = process_rank;
    cartesian_grid_communicator_ = cartesian_grid_communicator;
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      neighbours_ranks_[border] = neighbours_ranks[border];
      // Check array are valid for pointer arithmetics
      if (!borders_to_send(border).isStorageContiguous() ||
          !received_borders(border).isStorageContiguous()) {
        printf("Proc[%i]: Buffers for halo exchange is not contiguous!!!\n",
               process_rank_);
        return kErrorExchangeBufferIsNotContiguous;
      }  // end of if error: buffers are not contiguous
      if (borders_to_send(border).size() !=
          received_borders(border).size()) {
        printf("Proc[%i]: Buffers to send and recieve halo do not match!!!\n",
               process_rank_);
        return kErrorSendAndReceiveBuffersHasDifferentSizes;
      }  // end of if error: borders size send != receive
      borders_to_send_[border] = borders_to_send(border).data();
      received_borders_[border] = received_borders(border).data();
      number_of_elements_to_send_[border] = borders_to_send(border).size();
    }  // end of for border
    return kDone;
  }  // end of HaloToExchange::Init()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Check subdomains start and finish indexes to be
  /// correlated with neighbours indexes.
  ///
  /// Corespnding indexes of neighbours in selected dimension should
  /// have difference equal to one. Subomains at border of
  /// computational domain should have zero or domains length value.
  int HaloExchangeProcess::CheckSubdomainIndexes() {
    int64_t index_for_output[kDimensions*2];
    int64_t index_from_neighbour[kDimensions*2];
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      index_for_output[axis] = subdomain_start_index_[axis];
      index_for_output[axis+kDimensions] = subdomain_finish_index_[axis];
      index_from_neighbour[axis] = -1;
      index_from_neighbour[axis+kDimensions] = -1;
    }
    int send_size = 1;
    MPI_Request isend_request[kDimensions*2];
    MPI_Request irecv_request[kDimensions*2];
    MPI_Status status;
    /// @todo3 Check for MPI_INT64_T type, change to it MPI_LONG_LONG
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      int direction_tag = border < kDimensions
        ? border + kDimensions : border - kDimensions;
      MPI_Irecv(&index_from_neighbour[border], send_size, MPI_LONG_LONG,
                neighbours_ranks_[border], direction_tag,
                cartesian_grid_communicator_, &irecv_request[border]);
    }  // end of for borders receive
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      int direction_tag = border;
      MPI_Isend(&index_for_output[border], send_size, MPI_LONG_LONG,
                neighbours_ranks_[border], direction_tag,
                cartesian_grid_communicator_, &isend_request[border]);
    }  // end of for borders send
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (neighbours_ranks_[border] == MPI_PROC_NULL) continue;
      MPI_Wait(&isend_request[border], &status);
      MPI_Wait(&irecv_request[border], &status);
    }  // end of for borders wait
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      int64_t size = simulation_core_.simulation_input_config_.
        grid_input_config_.
        get_total_grid_length(static_cast<Axis>(border%kDimensions))-1;
      if (index_from_neighbour[border] != -1) {
        int64_t difference = abs(index_for_output[border]
                                 -index_from_neighbour[border]);
        if (difference !=1 && difference != size) {
          printf("Process[%i] border[%i]: Error! Wrong index difference!\n",
                 process_rank_, border);
          printf("Index %li from, %li for: != 1 or %li\n",
                 index_from_neighbour[border],
                 index_for_output[border],
                 size);
          return kErrorWrongIndexDifference;
        }  // end of if error
      } else {  // end of if border is domains internal
        if (border < kDimensions && index_for_output[border] != 0) {
          printf("Process[%i]: Error! Wrong start index!\n",
                 process_rank_);
          return kErrorWrongIndexDifference;
        }  // end of if nonzero index of most left (bottom, back) border
        if (border >= kDimensions && index_for_output[border] != size) {
          printf("Process[%i]: Error! Wrong finish index!\n",
                 process_rank_);
          return kErrorWrongIndexDifference;
        }  // end of if most right (top, front) index doesn`t
           // correlates with domain length
      }  // end of else sudmain`s border is domain external border
    }  // end of for borders check errors
    return kDone;
  }  // end of HaloExchangeProcess::CheckSubdomainIndexes()
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
    // Split for each axis independantly.
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      // Number of vertices.
      double total_size = simulation_core_.simulation_input_config_.
        grid_input_config_.get_total_grid_length(static_cast<Axis>(axis));
      // Number of equivalent vertices.
      double total_computational_size = total_size;
      // For axis with PML the computational size is bigger.
      // PML vertices are part of total grid.
      if (boundary_condition[axis] == kBoundaryConditionPML)
        total_computational_size += pml_width*(pml_computational_ratio - 1.0);
      if (boundary_condition[axis+kDimensions] == kBoundaryConditionPML)
        total_computational_size += pml_width*(pml_computational_ratio - 1.0);
      // @todo3 check if subdomain_computational_size is not too small (< 1).
      // Optimal number of equivalent vertices per subdomain dimension.
      double subdomain_computational_size =
        total_computational_size/static_cast<double>(subdomains_[axis]);
      // Evaluate grid coords from current MPI process coords.
      // Set undefined value.
      double subdomain_start_coord = -1;
      double subdomain_finish_coord = -1;
      // For wide PML with too many processes pure PML subdomains are
      // possible.
      double pml_subdomain_size = subdomain_computational_size
        /pml_computational_ratio;
      int64_t  pml_subdomains_number_from_left = 0;
      int64_t  pml_subdomains_number_from_right = 0;
      double transient_subdomain_pml_vertices_from_left = 0;
      double transient_subdomain_pml_vertices_from_right = 0;
      // Starting from first (left, bottom, back) side of domain.
      if (boundary_condition[axis] == kBoundaryConditionPML) {
        pml_subdomains_number_from_left
          = floor(pml_width/pml_subdomain_size);
        transient_subdomain_pml_vertices_from_left = pml_width -
          pml_subdomain_size * pml_subdomains_number_from_left;
      }  // end of if first boundary is PML
      // Second (right, top, front) side of domain.
      if (boundary_condition[axis+kDimensions] == kBoundaryConditionPML) {
        pml_subdomains_number_from_right
          = floor(pml_width/pml_subdomain_size);
        transient_subdomain_pml_vertices_from_right  = pml_width -
          pml_subdomain_size * pml_subdomains_number_from_right;
      }  // end of if second boundary is PML
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
        left_transient_subdomain_start_coord + subdomain_computational_size
        - transient_subdomain_pml_vertices_from_left
          * (pml_computational_ratio - 1.0);
      double right_transient_subdomain_start_coord =
        right_transient_subdomain_finish_coord - subdomain_computational_size
        + transient_subdomain_pml_vertices_from_right
          * (pml_computational_ratio - 1.0);
      // Special case: left and right PML are in same transient subdomain.
      if (left_transient_subdomain_finish_coord
          > right_transient_subdomain_finish_coord) {
        left_transient_subdomain_finish_coord
          = right_transient_subdomain_finish_coord;
        // printf("axis[%i] Special case : 1 \n", axis);
      }
      if (left_transient_subdomain_start_coord
          > right_transient_subdomain_start_coord) {
        right_transient_subdomain_start_coord =
          left_transient_subdomain_start_coord;
        // printf("axis[%i] Special case : 2 \n", axis);
      }
      if (my_coords_[axis] == pml_subdomains_number_from_left) {
        subdomain_start_coord = left_transient_subdomain_start_coord;
        subdomain_finish_coord = left_transient_subdomain_finish_coord;
      }  // end of if transient subdomain near first boundary
      if (subdomains_[axis] - 1 - my_coords_[axis]
          == pml_subdomains_number_from_right
          // Exclude special case: left and right PML are in same
          // transient subdomain.
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
      // HACK! try to deal with double rounding unstability.
      double ceil_subdomain_start_coord = ceil(subdomain_start_coord);
      if (abs(ceil_subdomain_start_coord - subdomain_start_coord) < 0.01
          && subdomain_start_coord > 0.1)
        subdomain_start_coord -= 0.1;
      double ceil_subdomain_finish_coord = ceil(subdomain_finish_coord);
      if (abs(ceil_subdomain_finish_coord - subdomain_finish_coord) < 0.01
          && abs(total_size - subdomain_finish_coord) > 0.1)
        subdomain_finish_coord -= 0.1;
      // Evaluate index!
      subdomain_start_index_[axis] = ceil(subdomain_start_coord);
      subdomain_finish_index_[axis] = ceil(subdomain_finish_coord) - 1;
      subdomain_size_[axis] = subdomain_finish_index_[axis]
        - subdomain_start_index_[axis] + 1;
      int halo_width = simulation_core_.simulation_input_config_.halo_width();
      if (subdomain_size_[axis] < halo_width && subdomain_size_[axis] != 1) {
        printf("Proc[%i] ERR! Subdomain size[%i] %li less than halo width %i\n",
               process_rank_, axis, subdomain_size_[axis], halo_width);
        printf("Use less MPI procs, smaller halo or bigger domain!\n");
        return kErrorSubdomainSizeLessThanHaloWidth;
      }
      /// @todo3 remove output in HaloExchangeProcess::EvaluateSubdomainSize().
      // int zeros = 0;
      // int non_zero_axis = -1;
      // for (int i = 0; i < kDimensions; ++i)
      //   if (my_coords_[i] == 0)
      //     zeros += 1;
      //   else
      //     non_zero_axis = i;
      // if ((zeros == 2 && axis == non_zero_axis) || zeros == 3) {
      //   printf("(%2i,%2i,%2i) axis[%i]: %3li -> %3li  (%3li): vert%7.2f\n",
      //          my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ],
      //          axis,
      //          subdomain_start_index_[axis],
      //          subdomain_finish_index_[axis],
      //          subdomain_size_[axis],
      //          subdomain_computational_size);
      // }  // end if process with coords output
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
      if (boundary_condition[axis] == kBoundaryConditionPML)
        length[axis] += ceil(pml_width*(pml_computational_ratio - 1.0));
      if (boundary_condition[axis+kDimensions] == kBoundaryConditionPML)
        length[axis] += ceil(pml_width*(pml_computational_ratio - 1.0));
    }  // end of for dimensions
    const int64_t all_processes = processes_total_number_;
    // (1) Topology decomposition.
    SimpleTopologyDecomposition3D(all_processes, length, subdomains_);
    if (process_rank_ == kOutput) {
      printf("Used MPI procs %li of %i (%.3g%%), subdomains_ (%li,%li,%li), ",
             MultiplyComponents(subdomains_), processes_total_number_,
             MultiplyComponents(subdomains_)
             /static_cast<double>(processes_total_number_)*100,
             subdomains_[kAxisX],
             subdomains_[kAxisY], subdomains_[kAxisZ]);
      printf("domain length (%li,%li,%li) equiv ->(%li,%li,%li)\n",
             orig_length[kAxisX],
             orig_length[kAxisY], orig_length[kAxisZ],
             length[kAxisX],
             length[kAxisY], length[kAxisZ]);
    }  // end of if process_rank_ == kOutput
    // (2) Start new communicator.
    if (StartCartesianGridCommunicator() != kDone)
      return kErrorProcessNotInGrid;
    // (3) Find neighbours.
    EvaluateCoordsAndRanks();
    /// @todo3 Remove output in HaloExchangeProcess::RunDecomposition()
    // printf("(%2i,%2i,%2i) rank: %2i  X(%2i,%2i) Y(%2i,%2i) Z(%2i,%2i)\n",
    //        my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ],
    //        process_rank_,
    //        neighbours_ranks_[kBorderLeft], neighbours_ranks_[kBorderRight],
    //        neighbours_ranks_[kBorderBottom], neighbours_ranks_[kBorderTop],
    //        neighbours_ranks_[kBorderBack],  neighbours_ranks_[kBorderFront]);
    // (4) Find current process subdomain size.
    EvaluateSubdomainSize();
    if (CheckSubdomainIndexes() != kDone) return kErrorWrongIndexDifference;
    return kDone;
  }  // end of HaloExchangeProcess::RunDecomposition
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Evaluate current process coords and neighbours ranks.
  int HaloExchangeProcess::EvaluateCoordsAndRanks() {
    MPI_Cart_coords(cartesian_grid_communicator_, process_rank_,
                    kDimensions, my_coords_);
    int upwards_shift = 1, downwards_shift = -1;
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      int rank_source, rank_dest;
      MPI_Cart_shift(cartesian_grid_communicator_, axis,
                     upwards_shift, &rank_source, &rank_dest);
      neighbours_ranks_[axis+kDimensions] = rank_dest;
      MPI_Cart_shift(cartesian_grid_communicator_, axis,
                     downwards_shift, &rank_source, &rank_dest);
      neighbours_ranks_[axis] = rank_dest;
    }  // end of for dimensions
    // Set periodical boundary conditions.
    int boundary_condition[kDimensions * 2];
    for (int border = kBorderLeft; border < kDimensions*2; ++border)
      boundary_condition[border] = simulation_core_.simulation_input_config_.
        boundary_condition(static_cast<BorderPosition>(border));
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      if (boundary_condition[border] != kBoundaryConditionPeriodical) continue;
      int axis = border % kDimensions;
      int direction_to_border = border < kDimensions ? -1 : 1;
      int max_shift = (subdomains_[axis] - 1);
      int rank_source, rank_dest;
      if (my_coords_[axis] == 0 && direction_to_border < 0) {
        MPI_Cart_shift(cartesian_grid_communicator_, axis,
                       max_shift, &rank_source, &rank_dest);
        neighbours_ranks_[border] = rank_dest;
        continue;
      }
      if (my_coords_[axis] == max_shift && direction_to_border > 0) {
        MPI_Cart_shift(cartesian_grid_communicator_, axis,
                       -max_shift, &rank_source, &rank_dest);
        neighbours_ranks_[border] = rank_dest;
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
  /// @breif Prepare simulation_core_ to start simulation.
  ///
  /// Initialize simulation_core_ member data.
  int HaloExchangeProcess::InitSimulation() {
    int init_status = simulation_core_.Init(subdomain_size_, process_rank_,
                                            neighbours_ranks_, my_coords_,
                                            subdomain_start_index_,
                                            subdomain_finish_index_);
    if (init_status != kDone) return init_status;
    simulation_core_.SetGridData();
    init_status = halo_to_exchange_.Init(simulation_core_.borders_to_send_,
                                         simulation_core_.received_borders_,
                                         process_rank_, neighbours_ranks_,
                                         cartesian_grid_communicator_);
    if (init_status != kDone) return init_status;
    return kDone;
  }  // end of HaloExchangeProcess::InitSimulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::Init(int argc, char *argv[]) {
    MPI_Comm_size(MPI_COMM_WORLD, &processes_total_number_);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    if (argc < 2) {
      printf("Proc[%i]: Error! Config file name was not defined!\n",
               process_rank_);
      return kErrorConfigFileNameWasNotDefined;
    }
    if (argc > 2) {
      printf("Proc[%i]: Error! Too many input parameters!\n",
               process_rank_);
      return kErrorConfigFileNameWasNotDefined;
    }
    std::string config_file_name(argv[1]);
    simulation_core_.simulation_input_config_.
      set_config_file_name(config_file_name);
    if (simulation_core_.simulation_input_config_.ReadConfig() != kDone)
      return kErrorTooWidePml;
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
    // Set floor and ceil normalized_length guess.
    for (int i = kAxisX; i < kDimensions; ++i) {
      // Reduced dimenstion is not optimized.
      if (is_dimension_reduced[i] == 1) {
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
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @breif Run simulation_core_ and halo exchange.
  ///
  /// This function runs FDTD simulation in MPI enviroment.
  int HaloExchangeProcess::RunSimulation() {
    if (simulation_core_.status() != kSimulationStatusInitiated) {
      printf("Process[%i] Error! Trying to run uninitiated simulation core\n",
             process_rank_);
      return kErrorUninitiatedSimulationCore;
    }  // end of if simulation core is not initiated
    double start_time, end_time;
    start_time = MPI_Wtime();
    while (simulation_core_.StepTime()) {
      int errors = 0;
      /// @todo1 change to for (a little bit faster)
      errors += timer_.Start("PrepareBordersToSend()");
      simulation_core_.PrepareBordersToSend();
      errors += timer_.Stop("PrepareBordersToSend()");
      // Sending borders, prepared with simulation_core. See also
      // HaloExchangeProcess::InitSimulation(). Borders are used
      // directly by pointer.
      errors += timer_.Start("StartNonBlockingExchange()");
      halo_to_exchange_.StartNonBlockingExchange();
      errors += timer_.Stop("StartNonBlockingExchange()");
      errors += timer_.Start("Snapshot()");
      simulation_core_.Snapshot();
      errors += timer_.Stop("Snapshot()");
      errors += timer_.Start("PrepareSource()");
      simulation_core_.PrepareSource();
      errors += timer_.Stop("PrepareSource()");
      errors += timer_.Start("FinishNonBlockingExchange()");
      halo_to_exchange_.FinishNonBlockingExchange();
      errors += timer_.Stop("FinishNonBlockingExchange()");
      errors += timer_.Start("DoBorderStep()");
      simulation_core_.DoBorderStep();
      errors += timer_.Stop("DoBorderStep()");
      errors += timer_.Start("DoStep()");
      /// @todo1 Replace DoStep() call with direct RunAlgorithm call
      simulation_core_.DoStep();
      errors += timer_.Stop("DoStep()");
      errors += timer_.Start("CycleSnapshots()");
      simulation_core_.CycleSnapshots();
      errors += timer_.Stop("CycleSnapshots()");
      if (errors) return kError;
      // Evaluate profile
    }  // end of while time is stepping
    end_time = MPI_Wtime();
    simulation_core_.Snapshot();
    double total_time_stepping = end_time-start_time;
    if (process_rank_ == kOutput) {
      printf("Stepping(%li) took %.2f seconds (%f s/step).\n",
             simulation_core_.total_time_steps(),
             total_time_stepping,
             total_time_stepping/simulation_core_.total_time_steps());
        timer_.PrintAllRelative(total_time_stepping);
      //  timer_.PrintAll();
      //  timer_.PrintByKey();
      double wasted = total_time_stepping - timer_.GetAllTotalTime();
      printf("*\t%04.1f%%\t%3.2f s\t%s\n",
             wasted/total_time_stepping*100.0, wasted,  "**wasted**");

    }
      // printf("DoStep() took %f seconds (%.1f%% from total).\n",
      //        do_step_total_time, do_step_total_time/total_time_stepping*100);
    return kDone;
  }  // end of HaloExchangeProcess::RunSimulation()
}  // end of namespace onza
