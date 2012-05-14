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
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::Init() {
    process_rank_ = MPI::COMM_WORLD.Get_rank();
    processes_total_number_ = MPI::COMM_WORLD.Get_size();
    int init_status = simulation_core_.Init();
    if (init_status != kDone) return init_status;
    return kDone;
  }  // end of HaloExchangeProcess::Init()
  /// @brief Run simulation domain decomposition
  ///
  /// With no information about cluster topology simply minimize surface
  /// of exchange border.
  ///
  ///  @todo Domain decomposition should be correlated with topology of cluster,
  ///  starting from intra node cores\processors layout.
  ///
  /// @return 0 or error code.
  int  HaloExchangeProcess::RunDecomposition() {
    if (simulation_core_.simulation_input_config_.
        status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;
    const int64_t length_x = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisX);
    const int64_t length_y = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisY);
    const int64_t length_z = simulation_core_.simulation_input_config_.
      grid_input_config_.get_total_grid_length(kAxisZ);
    // const int64_t all_processes = processes_total_number_;
    const int64_t all_processes = 9;  /// @todo Remove debug string.
    int64_t best_nx, best_ny, best_nz;
    StarTopologyDecomposition(all_processes, length_x, length_y, length_z,
                              best_nx, best_ny, best_nz);
    if (process_rank_ == 0) {
      printf("final all_processes = %li,  best_nx = %li, best_ny = %li, best_nz = %li\n",
             all_processes, best_nx, best_ny, best_nz);
    }  // end of if (process_rank_ == 0)
    return kDone;
  }  // end of HaloExchangeProcess::RunDecomposition
  /// @brief Optized decomposition for star topology of supercomputer.
  int HaloExchangeProcess::StarTopologyDecomposition(const int64_t all_processes,
                                                     const int64_t length_x,
                                                     const int64_t length_y,
                                                     const int64_t length_z,
                                                     int64_t &best_nx,
                                                     int64_t &best_ny,
                                                     int64_t &best_nz) {
    double single_cell_optimal_volume = static_cast<double>
      (length_x * length_y * length_z) / static_cast<double>(all_processes);
    double optimal_length = pow(single_cell_optimal_volume, 1.0/3.0);
    double nx = length_x/optimal_length;
    double ny = length_y/optimal_length;
    double nz = length_z/optimal_length;
    if (process_rank_ == 0) {
      printf("HaloExchange length_x = %li, length_x = %li, length_x = %li\n",
             length_x, length_y, length_z);
      printf("Total volume = %li, single cell volume = %.5g, cell length = %.5g\n",
             length_x * length_y * length_z, single_cell_optimal_volume, optimal_length);
      printf("all processes = %li,  nx = %g, ny = %g, nz = %g\n", all_processes, nx, ny, nz);
    }  // end of if (process_rank_ == 0)
    int64_t floor_nx = static_cast<int64_t>(floor(nx));
    if (floor_nx < 1) {floor_nx = 1;}
    if (floor_nx > all_processes) {floor_nx = all_processes;}
    int64_t ceil_nx = static_cast<int64_t>(ceil(nx));
    if (ceil_nx > all_processes) {ceil_nx = all_processes;}

    int64_t floor_ny = static_cast<int64_t>(floor(ny));
    if (floor_ny < 1) {floor_ny = 1;}
    if (floor_ny > all_processes) {floor_ny = all_processes;}
    int64_t ceil_ny = static_cast<int64_t>(ceil(ny));
    if (ceil_ny > all_processes) {ceil_ny = all_processes;}

    int64_t floor_nz = static_cast<int64_t>(floor(nz));
    if (floor_nz < 1) {floor_nz = 1;}
    if (floor_nz > all_processes) {floor_nz = all_processes;}
    int64_t ceil_nz = static_cast<int64_t>(ceil(nz));
    if (ceil_nz > all_processes) {ceil_nz = all_processes;}

    // const int variant_size = 4;
    // int64_t nx_variant[] = {floor_nx, ceil_nx,
    //                         ceil_nx+1 > all_processes ? all_processes : ceil_nx+1,
    //                         ceil_nx+2 > all_processes ? all_processes : ceil_nx+2};
    // int64_t ny_variant[] = {floor_ny, ceil_ny,
    //                         ceil_ny+1 > all_processes ? all_processes : ceil_ny+1,
    //                         ceil_ny+2 > all_processes ? all_processes : ceil_ny+2};
    // int64_t nz_variant[] = {floor_nz, ceil_nz,
    //                         ceil_nz+1 > all_processes ? all_processes : ceil_nz+1,
    //                         ceil_nz+2 > all_processes ? all_processes : ceil_nz+2};
    const int variant_size = 5;
    int64_t nx_variant[] = {floor_nx-1 < 1 ? 1 : floor_nx-1,
                            floor_nx, ceil_nx,
                            ceil_nx+1 > all_processes ? all_processes : ceil_nx+1,
                            ceil_nx+2 > all_processes ? all_processes : ceil_nx+2};
    int64_t ny_variant[] = {floor_ny-1 < 1 ? 1 : floor_ny-1,
                            floor_ny, ceil_ny,
                            ceil_ny+1 > all_processes ? all_processes : ceil_ny+1,
                            ceil_ny+2 > all_processes ? all_processes : ceil_ny+2};
    int64_t nz_variant[] = {floor_nz-1 < 1 ? 1 : floor_nz-1,
                            floor_nz, ceil_nz,
                            ceil_nz+1 > all_processes ? all_processes : ceil_nz+1,
                            ceil_nz+2 > all_processes ? all_processes : ceil_nz+2};
    best_nx = nx_variant[0];
    best_ny = ny_variant[0];
    best_nz = nz_variant[0];
    int64_t best_size = best_nx * best_ny * best_nz;
    int64_t best_norm = (best_nx-1)*length_y*length_z
                      + (best_ny-1)*length_x*length_z
                      + (best_nz-1)*length_x*length_y;
    for (int i = 0; i < variant_size; ++i) {
      for (int j = 0; j < variant_size; ++j) {
        for (int k = 0; k < variant_size; ++k) {
          int64_t new_size = nx_variant[i] * ny_variant[j] * nz_variant[k];
          int64_t new_norm = (nx_variant[i]-1)*length_y*length_z
                           + (ny_variant[j]-1)*length_x*length_z
                           + (nz_variant[k]-1)*length_x*length_y;
//           if (process_rank_ == 0) {
//             printf("testing new_size = %li, best_size = %li, new norm = %li, best_norm = %li, \
// nx_variant = %li, ny_variant = %li, nz_variant = %li\n",
//                    new_size, best_size, new_norm, best_norm,
//                    nx_variant[i], ny_variant[j], nz_variant[k]);
//           }
          if (new_size > all_processes) continue;
          if (new_size == best_size && new_norm < best_norm) {
            best_norm = new_norm;
            best_nx = nx_variant[i];
            best_ny = ny_variant[j];
            best_nz = nz_variant[k];
          }
          if (new_size > best_size) {
            best_size = new_size;
            best_norm = new_norm;
            best_nx = nx_variant[i];
            best_ny = ny_variant[j];
            best_nz = nz_variant[k];
          }  // end of new best values
        }
      }  // end of for j
    }  // end of for i
    if (process_rank_ == 0) {
      printf("pre final all_processes = %li,  best_nx = %li, best_ny = %li, best_nz = %li\n",
             all_processes, best_nx, best_ny, best_nz);
    }  // end of if (process_rank_ == 0)
    return kDone;
  }  // end of HaloExchangeProcess::StarTopologyDecomposition

}  // end of namespace onza

