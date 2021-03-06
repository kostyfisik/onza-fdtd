///
/// @file   basic-fdtd.cc
/// @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
/// @copyright 2012 Konstantin Ladutenko
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
/// @date   Wed May  2 13:42:32 2012
///
/// @brief  Simulation data and stepping algorithm
#include <mpi.h>
#include <blitz/array.h>
//debug for blitz::Array output
// #include <iostream>
#include <map>
#include <cstdio>
#include "./basic-fdtd.h"
#include "../common.h"
namespace onza {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Prepare borders to send
  ///
  /// Copy grid data_ slices to borders_to_send_
  void BasicSimulationCore::PrepareBordersToSend() {
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      for (int component = 0;
           component < number_of_components_to_exchange_;
           ++component) {
        // Get slice of data_ for current border and current
        // component.  In fist (component) dimension it will have the
        // only element (current component). Copy slice (reduced to 3D
        // array by component index) from data_ to corresponding slice of
        // borders_to_send_.
        const int the_only_component_index = 0;
        int opposite_border = (border + kDimensions) % (kDimensions*2);
        if (neighbours_ranks_[border] != MPI_PROC_NULL) {
          if (neighbours_ranks_[border] != process_rank_)
            borders_to_send_(border)(component, all_x_, all_y_, all_z_)
              = data_(borders_to_send_range_(border, component))
              (the_only_component_index, all_x_, all_y_, all_z_);
          else
            received_borders_(opposite_border)
              (component, all_x_, all_y_, all_z_)
              = data_(borders_to_send_range_(border, component))
              (the_only_component_index, all_x_, all_y_, all_z_);
        }  // end of if neighbour receiver if valid
      }  // end of for component
    }  // end of for border
  }  // end of BasicSimulationCore::PrepareBordersToSend()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do time stepping for current subdomain.
  ///
  /// @warning  It is repeated MANY times!
  int BasicSimulationCore::StepTime() {
    // Step with data_ is done. Stepping  time...
    ++local_time_step_;
    if (total_time_steps_ < local_time_step_) return  kSimulationStatusFinished;
    return kSimulationStatusRunning;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Snapshot component from data_ to field.
  ///
  /// @warning  It is repeated MANY times!
  void BasicSimulationCore::Snapshot() {
    if (snapshot_interval_ * snapshot_frame_ < local_time_step_) {
      ++snapshot_frame_;
      char filename[100];
      snprintf(filename, sizeof(filename),
               "Ez(x)[%i]-time%.7li-domain(x%.3dy%.3dz%.3d).onza",
               algorithm_, local_time_step_,
               my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ]);
      FILE *snapshot;
      snapshot = fopen(filename, "w");
      // int offset_x = my_coords_[kAxisX] * subdomain_size_[kAxisX];
      /// @todo2 make snapshot output  call be warper of algorithm
      /// dependent function call (call snapshot  with function pointer
      /// assigned at init time)
      switch (algorithm_) {
      case kAlgorithmSimpleX1D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot, "%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  data_(static_cast<int>(kEz), i, 0, 0));
        break;
      case kAlgorithmSimpleY1D:
        for (int i = 0; i < subdomain_size_[kAxisY]; ++i)
          fprintf(snapshot, "%li\t%g\n", i+subdomain_start_index_[kAxisY],
                  data_(static_cast<int>(kEz), 0, i, 0));
        break;
      case kAlgorithmSimpleZ1D:
        for (int i = 0; i < subdomain_size_[kAxisZ]; ++i)
          fprintf(snapshot, "%li\t%g\n", i+subdomain_start_index_[kAxisZ],
                  data_(static_cast<int>(kEz), 0, 0, i));
        break;
      case kAlgorithmSimpleTMz2D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot, "%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  data_(static_cast<int>(kEz), i,
                        static_cast<int>(subdomain_size_[kAxisY]/2),
                        static_cast<int>(subdomain_size_[kAxisZ]/2)));
        break;
      case kAlgorithmSimple3D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot, "%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  data_(static_cast<int>(kEz), i,
                        static_cast<int>(subdomain_size_[kAxisY]/2),
                        static_cast<int>(subdomain_size_[kAxisZ]/2)));
        break;
        // default:  // @todo3 After moving snapshot call to warper
                     // we can handle error there
        //   printf("Error! Should use some FDTD algorithm!\n");
        //   return kErrorWrongAlgorithm;
      }  // end of switch algorithm
      fclose(snapshot);
    }  // end of if it is time to output snapshot
  }  // end of BasicSimulationCore::Snapshot()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Prepare source component in data_.
  ///
  /// @warning  It is repeated MANY times!
  void BasicSimulationCore::PrepareSource() {
    //debug 3D
    /// @todo2 make source preparation call be warper of algorithm
    /// dependent function call (call prepare with function pointer
    /// assigned at init time)
    //debug
    if (algorithm_ == kAlgorithmSimple3D) {
      data_(static_cast<int>(kSrcEz), all_x_, all_y_, all_z_) =
        exp(-pow2(local_time_step_ - 50.) / 200.);
        return;
    }  // end if special benchmarck source in each MPI process
    //debug
    if (algorithm_ == kAlgorithmSimpleTMz2D) {
      data_(static_cast<int>(kSrcEz), all_x_, all_y_, all_z_) =
        exp(-pow2(local_time_step_ - 50.) / 200.);
        return;
    }  // end if special benchmarck source in each MPI process
    //debug
    int64_t source_global_x = 50, source_global_y = 0, source_global_z = 0;
    if (subdomain_start_index_[kAxisX] <= source_global_x
        && subdomain_finish_index_[kAxisX] >= source_global_x
        && subdomain_start_index_[kAxisY] <= source_global_y
        && subdomain_finish_index_[kAxisY] >= source_global_y
        && subdomain_start_index_[kAxisZ] <= source_global_z
        && subdomain_finish_index_[kAxisZ] >= source_global_z) {
    // if (process_rank_ == 0) {
      int source_local_x = source_global_x - subdomain_start_index_[kAxisX];
      int source_local_y = source_global_y - subdomain_start_index_[kAxisY];
      int source_local_z = source_global_z - subdomain_start_index_[kAxisZ];
      switch (algorithm_) {
      case kAlgorithmSimpleX1D:
        data_(kSrcEz, source_local_x, all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleY1D:
        data_(kSrcEz, all_x_, source_local_y, all_z_) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleZ1D:
        data_(kSrcEz, all_x_, all_y_, source_local_z) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleTMz2D:
        //debug
        // data_(static_cast<int>(kSrcEz), 110, 140, all_z_) =
        //   exp(-pow2(local_time_step_ - 50.) / 200.);
        data_(static_cast<int>(kSrcEz), source_local_x, all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 50.) / 200.);
        break;
      case kAlgorithmSimple3D:
        //debug point
        // data_(static_cast<int>(kSrcEz),source_local_x, source_local_y,
        //       source_local_z) = exp(-pow2(local_time_step_ - 300.) / 4000.);
        //debug line
        // data_(static_cast<int>(kSrcEz),
        //       source_local_x, source_local_y, all_z_) =
        //   exp(-pow2(local_time_step_ - 100.) / 400.);
        // plane
        data_(static_cast<int>(kSrcEz),
              static_cast<int>(source_local_x), all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 100.) / 400.);
        break;
      // default:
      //   printf("Error! Should use some FDTD algorithm!\n");
      //   return kErrorWrongAlgorithm;
      }  // end of switch algorithm
    }  //  end of if source is in subdomain
  }  // end of BasicSimulationCore::PrepareSource()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for 3D case.
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::AlgorithmSimple3D(blitz::Range x,
                                           blitz::Range y,
                                           blitz::Range z) {
    // current time snapshot in snapshots array
    int t = time_depth_ - 2;
    // constant shift to normalize equation to space range.
    int c = 0;  /// @todo1 use const int c0, c1, cm1 , c2 and so on.
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    // Update H
    // Hx
    data_snapshot_        (1+t)(kHx  , x    , y    , z) =
      data_snapshot_        (t)(kChxh, x    , y    , z)
      * data_snapshot_      (t)(kHx  , x    , y    , z)
      +
      data_snapshot_        (t)(kChxe, x    , y    , z)
      * ((data_snapshot_    (t)(kEy  , x    , y    , z + 1)
          - data_snapshot_  (t)(kEy  , x    , y    , z))
         -
         (data_snapshot_    (t)(kEz  , x    , y + 1, z)
          - data_snapshot_  (t)(kEz  , x    , y    , z)));
    // Hy
    data_snapshot_        (1+t)(kHy  , x    , y    , z) =
      data_snapshot_        (t)(kChyh, x    , y    , z)
      * data_snapshot_      (t)(kHy  , x    , y    , z)
      +
      data_snapshot_        (t)(kChye, x    , y    , z)
      * ((data_snapshot_    (t)(kEz  , x + 1, y    , z)
          - data_snapshot_  (t)(kEz  , x    , y    , z))
         -
         (data_snapshot_    (t)(kEx  , x    , y    , z + 1)
          - data_snapshot_  (t)(kEx  , x    , y    , z)));
    // Hz
    data_snapshot_        (1+t)(kHz  , x    , y    , z) =
      data_snapshot_        (t)(kChzh, x    , y    , z)
      * data_snapshot_      (t)(kHz  , x    , y    , z)
      +
      data_snapshot_        (t)(kChze, x    , y    , z)
      * ((data_snapshot_    (t)(kEx  , x    , y + 1, z)
          - data_snapshot_  (t)(kEx  , x    , y    , z))
         -
         (data_snapshot_    (t)(kEy  , x + 1, y    , z)
          - data_snapshot_  (t)(kEy  , x    , y    , z)));
    // Update E
    // Ex
    c = 1;
    data_snapshot_        (1+t)(kEx  , x    , y     + c, z     + c) =
      data_snapshot_        (t)(kCexe, x    , y     + c, z     + c)
      * data_snapshot_      (t)(kEx  , x    , y     + c, z     + c)
      +
      data_snapshot_        (t)(kCexh, x    , y     + c, z     + c)
      * ((data_snapshot_  (1+t)(kHz  , x    , y     + c, z     + c)
          - data_snapshot_(1+t)(kHz  , x    , y - 1 + c, z     + c))
         -
         (data_snapshot_  (1+t)(kHy  , x    , y     + c, z     + c)
          - data_snapshot_(1+t)(kHy  , x    , y     + c, z - 1 + c)));
    // Ey
    data_snapshot_        (1+t)(kEy  , x     + c, y    , z     + c) =
      data_snapshot_        (t)(kCeye, x     + c, y    , z     + c)
      * data_snapshot_      (t)(kEy  , x     + c, y    , z     + c)
      +
      data_snapshot_        (t)(kCeyh, x     + c, y    , z     + c)
      * ((data_snapshot_  (1+t)(kHx  , x     + c, y    , z     + c)
          - data_snapshot_(1+t)(kHx  , x     + c, y    , z - 1 + c))
         -
         (data_snapshot_  (1+t)(kHz  , x     + c, y    , z     + c)
          - data_snapshot_(1+t)(kHz  , x - 1 + c, y    , z     + c)));
    // Ez
    data_snapshot_        (1+t)(kEz   , x     + c, y     + c, z) =
      data_snapshot_        (t)(kCeze , x     + c, y     + c, z)
      * data_snapshot_      (t)(kEz   , x     + c, y     + c, z)
      +
      data_snapshot_        (t)(kCezh , x     + c, y     + c, z)
      * ((data_snapshot_  (1+t)(kHy   , x     + c, y     + c, z)
          - data_snapshot_(1+t)(kHy   , x - 1 + c, y     + c, z))
         -
         (data_snapshot_  (1+t)(kHx   , x     + c, y     + c, z)
          - data_snapshot_(1+t)(kHx   , x     + c, y - 1 + c, z)))
      +
      data_snapshot_        (t)(kSrcEz, x     + c, y     + c, z);
  }  // end of BasicSimulationCore::AlgorithmSimple3D()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for TMz 2D case.
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::AlgorithmSimpleTMz2D(blitz::Range x,
                                                 blitz::Range y,
                                                 blitz::Range z) {
    // current time snapshot in snapshots array
    int t = time_depth_ - 2;
    // constant shift to normalize equation to space range.
    int c = 0;  /// @todo1 use const int c0, c1, cm1 , c2 and so on.
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    // Update H
    // Hx
    data_snapshot_        (1+t)(kHx  , x    , y    , z) =
      data_snapshot_        (t)(kChxh, x    , y    , z)
      * data_snapshot_      (t)(kHx  , x    , y    , z)
      -
      data_snapshot_        (t)(kChxe, x    , y    , z)
      * (data_snapshot_     (t)(kEz  , x    , y + 1, z)
         - data_snapshot_   (t)(kEz  , x    , y    , z));
    // Hy
    data_snapshot_        (1+t)(kHy  , x    , y    , z) =
      data_snapshot_        (t)(kChyh, x    , y    , z)
      * data_snapshot_      (t)(kHy  , x    , y    , z)
      +
      data_snapshot_        (t)(kChye, x    , y    , z)
      * (data_snapshot_     (t)(kEz  , x + 1, y    , z)
         - data_snapshot_   (t)(kEz  , x    , y    , z));
    // Update E
    // Ez
    c = 1;
    data_snapshot_        (1+t)(kEz   , x     + c, y     + c, z) =
      data_snapshot_        (t)(kCeze , x     + c, y     + c, z)
      * data_snapshot_      (t)(kEz   , x     + c, y     + c, z)
      +
      data_snapshot_        (t)(kCezh , x     + c, y     + c, z)
      * ((data_snapshot_  (1+t)(kHy   , x     + c, y     + c, z)
          - data_snapshot_(1+t)(kHy   , x - 1 + c, y     + c, z))
         -
         (data_snapshot_  (1+t)(kHx   , x     + c, y     + c, z)
          - data_snapshot_(1+t)(kHx   , x     + c, y - 1 + c, z)))
      +
      data_snapshot_        (t)(kSrcEz, x     + c, y     + c, z);
  }  // end of BasicSimulationCore::AlgorithmSimpleTMz2D()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for 1D case (length in kAxisX dimenstion).
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::AlgorithmSimpleX1D(blitz::Range x,
                                           blitz::Range y,
                                           blitz::Range z) {
    /// @todo1 Implement variable sized array for constants.
    // impedance of free space.
    double imp0 = 377.0;
    double inv_imp0 = 1.0/imp0;
    // current time snapshot in snapshots array
    int t = time_depth_ - 2;
    // constant shift to normalize equation to space range.
    int c = 0;   /// @todo1 use const int c0, c1, cm1 , c2 and so on.
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (1+t)(kHy, x     + c, y, z) =
      data_snapshot_    (t)(kHy, x     + c, y, z)
      + data_snapshot_  (t)(kEz, x + 1 + c, y, z) * inv_imp0
      - data_snapshot_  (t)(kEz, x     + c, y, z) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (1+t)(kEz    , x     + c, y, z) =
        data_snapshot_    (t)(kEz    , x     + c, y, z)
      + data_snapshot_  (1+t)(kHy    , x     + c, y, z)
        * data_snapshot_  (t)(kCezh, x     + c, y, z)
        * imp0
      - data_snapshot_  (1+t)(kHy    , x - 1 + c, y, z)
        * data_snapshot_  (t)(kCezh, x     + c, y, z)
        * imp0
      + data_snapshot_    (t)(kSrcEz , x     + c, y, z);
  }  // end of BasicSimulationCore::AlgorithmSimpleX1D()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for 1D case (length in kAxisY dimenstion).
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::AlgorithmSimpleY1D(blitz::Range x,
                                           blitz::Range y,
                                           blitz::Range z) {
    /// @todo1 Implement variable sized array for constants.
    // impedance of free space.
    double imp0 = 377.0;
    double inv_imp0 = 1.0/imp0;
    // current time snapshot in snapshots array
    int t = time_depth_ - 2;
    // constant shift to normalize equation to space range.
    int c = 0;  /// @todo1 use const int c0, c1, cm1 , c2 and so on.
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (1+t)(kHy, x, y     + c, z) =
      data_snapshot_    (t)(kHy, x, y     + c, z)
      + data_snapshot_  (t)(kEz, x, y + 1 + c, z) * inv_imp0
      - data_snapshot_  (t)(kEz, x, y     + c, z) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (1+t)(kEz    , x, y     + c, z) =
        data_snapshot_    (t)(kEz    , x, y     + c, z)
      + data_snapshot_  (1+t)(kHy    , x, y     + c, z)
        * data_snapshot_  (t)(kCezh, x, y     + c, z)
        * imp0
      - data_snapshot_  (1+t)(kHy    , x, y - 1 + c, z)
        * data_snapshot_  (t)(kCezh, x, y     + c, z)
        * imp0
      + data_snapshot_    (t)(kSrcEz , x, y     + c, z);
  }  // end of BasicSimulationCore::AlgorithmSimpleY1D()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for 1D case (length in kAxisZ dimenstion).
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::AlgorithmSimpleZ1D(blitz::Range x,
                                           blitz::Range y,
                                           blitz::Range z) {
    /// @todo1 Implement variable sized array for constants.
    // impedance of free space.
    double imp0 = 377.0;
    double inv_imp0 = 1.0/imp0;
    // current time snapshot in snapshots array
    int t = time_depth_ - 2;
    // constant shift to normalize equation to space range.
    int c = 0;  /// @todo1 use const int c0, c1, cm1 , c2 and so on.
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (1+t)(kHy, x, y, z     + c) =
      data_snapshot_    (t)(kHy, x, y, z     + c)
      + data_snapshot_  (t)(kEz, x, y, z + 1 + c) * inv_imp0
      - data_snapshot_  (t)(kEz, x, y, z     + c) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (1+t)(kEz    , x, y, z     + c) =
        data_snapshot_    (t)(kEz    , x, y, z     + c)
      + data_snapshot_  (1+t)(kHy    , x, y, z     + c)
        * data_snapshot_  (t)(kCezh, x, y, z     + c)
        * imp0
      - data_snapshot_  (1+t)(kHy    , x, y, z - 1 + c)
        * data_snapshot_  (t)(kCezh, x, y, z     + c)
        * imp0
      + data_snapshot_    (t)(kSrcEz , x, y, z     + c);
  }  // end of BasicSimulationCore::AlgorithmSimpleZ1D()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Cycle snapshots before new timestep.
  ///
  /// Rearrange references to snapshot array so that
  /// data_snapshot_(timestep) has data of data_snapshot_(timestep+1)
  /// @warning It is repeated MANY times!
  void BasicSimulationCore::CycleSnapshots() {
    data_.reference(data_snapshot_(0));
    for (int time = 0; time < time_depth_-1; ++time)
      data_snapshot_(time).reference(data_snapshot_(time+1));
    data_snapshot_(time_depth_-1).reference(data_);
    data_.reference(data_snapshot_(time_depth_-2));
  }  // end of BasicSimulationCore::CycleSnapshots()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for border part of grid data.
  ///
  /// @see PrepareBordersToSend().
  /// @warning  It is repeated MANY times!
  void BasicSimulationCore::DoBorderStep() {
    // Copy received border do grid data.
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      for (int component = 0;
           component < number_of_components_to_exchange_;
           ++component) {
        const int the_only_component_index = 0;
        // int opposite_border = (border + kDimensions) % (kDimensions*2);
        if (neighbours_ranks_[border] != MPI_PROC_NULL)
          data_(received_borders_range_(border, component))
                  (the_only_component_index, all_x_, all_y_, all_z_)
            = received_borders_(border)(component, all_x_, all_y_, all_z_);
      }  // end of for component
    }  // end of for border
    //debug
    // if (neighbours_ranks_[kBorderLeft] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(data_border_range_[kBorderLeft], inner_y_,
    //                         inner_z_);
    // if (neighbours_ranks_[kBorderRight] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(data_border_range_[kBorderRight], inner_y_,
    //                         inner_z_);
    // if (neighbours_ranks_[kBorderBottom] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(inner_x_, data_border_range_[kBorderBottom],
    //                         inner_z_);
    // if (neighbours_ranks_[kBorderTop] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(inner_x_, data_border_range_[kBorderTop],
    //                         inner_z_);
    // if (neighbours_ranks_[kBorderBack] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(inner_x_, inner_y_,
    //                         data_border_range_[kBorderBack]);
    // if (neighbours_ranks_[kBorderFront] != MPI_PROC_NULL)
    //   (this->*RunAlgorithm)(inner_x_, inner_y_,
    //                         data_border_range_[kBorderFront]);
  }  // end of BasicSimulationCore::DoStep()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for internal part of grid data.
  ///
  /// @warning  It is repeated MANY times!
  void BasicSimulationCore::DoStep() {
    (this->*RunAlgorithm)(inner_x_, inner_y_, inner_z_);
  }  // end of BasicSimulationCore::DoStep()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize grid data for all components
  int BasicSimulationCore::SetGridData() {
    // set all fields (Ex, Ey, Ez, Hx, Hy, Hz) to be zero. And everything else.
    data_ = 0;
    double epsilon = 1.0;
    // double epsilon = 2.1;
    /// @todo1 Implement variable sized array for constants.
    // impedance of free space.
    double imp0 = 377.0;
    double courant_3D = 1.0/sqrt(3.0);
    double courant_2D = 1.0/sqrt(2.0);
    switch (algorithm_) {
    case kAlgorithmSimpleX1D:
    case kAlgorithmSimpleY1D:
    case kAlgorithmSimpleZ1D:
      //debug
      data_(kCezh, all_x_, all_y_, all_z_) = 1.0/epsilon;
      break;
    case kAlgorithmSimpleTMz2D:
      //debug
      data_(kCeze, all_x_, all_y_, all_z_) = 1.0;
      data_(kCezh, all_x_, all_y_, all_z_) = courant_2D * imp0;
      data_(kChxh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChxe, all_x_, all_y_, all_z_) = courant_2D / imp0;
      data_(kChyh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChye, all_x_, all_y_, all_z_) = courant_2D / imp0;
      break;
    case kAlgorithmSimple3D:
      //debug
      data_(kCexe, all_x_, all_y_, all_z_) = 1.0;
      data_(kCexh, all_x_, all_y_, all_z_) = courant_3D * imp0;
      data_(kCeye, all_x_, all_y_, all_z_) = 1.0;
      data_(kCeyh, all_x_, all_y_, all_z_) = courant_3D * imp0;
      data_(kCeze, all_x_, all_y_, all_z_) = 1.0;
      data_(kCezh, all_x_, all_y_, all_z_) = courant_3D * imp0;
      data_(kChxh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChxe, all_x_, all_y_, all_z_) = courant_3D / imp0;
      data_(kChyh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChye, all_x_, all_y_, all_z_) = courant_3D / imp0;
      data_(kChzh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChze, all_x_, all_y_, all_z_) = courant_3D / imp0;
      break;
    default:
      throw std::invalid_argument("Error! Should use some FDTD algorithm!\n");
    }
    // Prepare snapshots. data_snapshot_(time_depth_ - 2) and data_
    // references to the same memory area.
    for (int i = 0; i < time_depth_; ++i) {
      if (i == time_depth_ - 2) continue;
      data_snapshot_(i) = data_snapshot_(time_depth_-2);
    }
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  BasicSimulationCore::Init(const int64_t subdomain_size[],
                                 int process_rank, int neighbours_ranks[],
                                 int my_coords[],
                                 int64_t subdomain_start_index[],
                                 int64_t subdomain_finish_index[]) {
    if (simulation_input_config_.status() != kInputConfigAllDone)
      throw std::invalid_argument("Error! Using input config too early!");
    status_ = kSimulationStatusUndefined;
    // Init member in alphabetical order
    algorithm_ = simulation_input_config_.algorithm();
    switch (algorithm_) {
    case kAlgorithmSimpleX1D:
      RunAlgorithm = &BasicSimulationCore::AlgorithmSimpleX1D;
      break;
    case kAlgorithmSimpleY1D:
      RunAlgorithm = &BasicSimulationCore::AlgorithmSimpleY1D;
      break;
    case kAlgorithmSimpleZ1D:
      RunAlgorithm = &BasicSimulationCore::AlgorithmSimpleZ1D;
      break;
    case kAlgorithmSimpleTMz2D:
      RunAlgorithm = &BasicSimulationCore::AlgorithmSimpleTMz2D;
      break;
    case kAlgorithmSimple3D:
      RunAlgorithm = &BasicSimulationCore::AlgorithmSimple3D;
      break;
    default:
      throw std::invalid_argument("Error! Should use some FDTD algorithm!\n");
    }  // end of switch algorithm
    halo_width_ = simulation_input_config_.halo_width();
    if (halo_width_ < 1) {
      throw std::invalid_argument("Error! FDTD algorithm`s should have some halo!\n");
    }  // end of if error
    local_time_step_ = 0;
    for (int border = kBorderLeft; border < kDimensions*2; ++border)
      neighbours_ranks_[border] = neighbours_ranks[border];
    number_of_components_to_exchange_ = simulation_input_config_.
      number_of_components_to_exchange();
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = simulation_input_config_.components_to_exchange();
    number_of_grid_data_components_ = simulation_input_config_.
      number_of_grid_data_components();
    process_rank_ = process_rank;
    snapshot_interval_ = simulation_input_config_.snapshot_interval();
    snapshot_frame_ = 0;
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      my_coords_[axis] = my_coords[axis];
      subdomain_size_[axis] = subdomain_size[axis];
      subdomain_start_index_[axis] = subdomain_start_index[axis];
      subdomain_finish_index_[axis] = subdomain_finish_index[axis];
    }  // end of for axis init
    time_depth_ = simulation_input_config_.time_depth();
    total_time_steps_ = simulation_input_config_.total_time_steps();
    // Resize grid data (in space and time) to fit subdomain of
    // current MPI process.
    int max_x = subdomain_size_[kAxisX]-1;
    int max_y = subdomain_size_[kAxisY]-1;
    int max_z = subdomain_size_[kAxisZ]-1;
    if (time_depth_ < 2) {
      throw std::invalid_argument("Error! FDTD algorithm`s time depth should be >= 2!\n");
    }  // end of if error
    if (!max_x && !max_y && !max_z) {
      printf("\tUse bigger model or less processes!\n");
      throw std::invalid_argument("Error! Single element subdomains are not supported!\n");
    }  // end of if error

    data_snapshot_.resize(time_depth_);
    // Mark reduce dimensions.
    int is_reduced_x = max_x == 0 ? 0 : 1;
    int is_reduced_y = max_y == 0 ? 0 : 1;
    int is_reduced_z = max_z == 0 ? 0 : 1;
    int halo_width_x = halo_width_ * is_reduced_x;
    int halo_width_y = halo_width_ * is_reduced_y;
    int halo_width_z = halo_width_ * is_reduced_z;
    for (int time = 0; time < time_depth_; ++time) {
      // Reduced dimensions are excluded from halo exchange to allow
      // pure 1D and 2D simulations with minimal overhead. While
      // reading config boundary conditions kBoundaryConditionReduced
      // should be set in this dimension.
      data_snapshot_(time).resize(blitz::shape(number_of_grid_data_components_,
               max_x + 1 + halo_width_x * 2,   max_y + 1 + halo_width_y * 2,
               max_z + 1 + halo_width_z * 2));
      data_snapshot_(time).reindexSelf(blitz::shape(0, -halo_width_x,
                                 -halo_width_y, -halo_width_z));
    }  // end of for time of snapshot
    data_.resize(blitz::shape(number_of_grid_data_components_,
                              max_x + 1 + halo_width_x * 2,
                              max_y + 1 + halo_width_y * 2,
                              max_z + 1 + halo_width_z * 2));
    data_.reindexSelf(blitz::shape(0, -halo_width_x,
                                   -halo_width_y, -halo_width_z));
    // Ranges for all data spacial coords.
    all_x_ = blitz::Range::all();
    all_y_ = blitz::Range::all();
    all_z_ = blitz::Range::all();
    // Special reference to current time step data.
    data_.reference(data_snapshot_(time_depth_-2));
    // Specify ranges in each dimension for FDTD time stepping. For
    // ease of defining FDTD equations range can be shifted in forward
    // or backwrd direction for any decimal value less or equal
    // halo_width_. Note that aligment shift for each FDTD equation
    // should be applied to provide all ranges to be valid.
    //debug
    inner_x_ = blitz::Range(-halo_width_x, max_x);
    inner_y_ = blitz::Range(-halo_width_y, max_y);
    inner_z_ = blitz::Range(-halo_width_z, max_z);
    // inner_x_ = blitz::Range(0, max_x);
    // inner_y_ = blitz::Range(0, max_y);
    // inner_z_ = blitz::Range(0, max_z);
    // Specify data range for each border. Due to choice of
    // aligment shift in FDTD equation low and hi index ranges have
    // different number of elements.
    data_border_range_[kBorderLeft] = blitz::Range(-halo_width_x,
                                                   halo_width_x);
    data_border_range_[kBorderBottom] = blitz::Range(-halo_width_y,
                                                     halo_width_y);
    data_border_range_[kBorderBack] = blitz::Range(-halo_width_z,
                                                   halo_width_z);
    data_border_range_[kBorderRight] = blitz::Range(max_x - halo_width_x,
                                                    max_x);
    data_border_range_[kBorderTop] = blitz::Range(max_y - halo_width_y,
                                                  max_y);
    data_border_range_[kBorderFront] = blitz::Range(max_z - halo_width_z,
                                                    max_z);
    // Specify borders ranges for halo exchange.
    // Used for quick selection in BasicSimulationCore::PrepareBordersToSend().
    borders_to_send_range_.resize(kDimensions*2,
                                  number_of_components_to_exchange_);
    received_borders_range_.resize(kDimensions*2,
                                   number_of_components_to_exchange_);
    for (int component = 0;
         component < number_of_components_to_exchange_; ++component) {
      // Try to shorten notation
      typedef blitz::RectDomain<1+kDimensions> rd;
      typedef blitz::TinyVector<int, 4> vec;
      int c = components_to_exchange_(component);
      // Ranges for border to be send.
      borders_to_send_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c, 0            , 0    , 0),
             vec(c, halo_width_-1, max_y, max_z));
      borders_to_send_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c, max_x - (halo_width_-1), 0    , 0),
             vec(c, max_x                  , max_y, max_z));
      borders_to_send_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c, 0    , 0            , 0),
             vec(c, max_x, halo_width_-1, max_z));
      borders_to_send_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c, 0    , max_y - (halo_width_-1), 0),
             vec(c, max_x, max_y                  , max_z));
      borders_to_send_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c, 0    , 0    , 0),
             vec(c, max_x, max_y, halo_width_-1));
      borders_to_send_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c, 0    , 0    , max_z - (halo_width_-1)),
             vec(c, max_x, max_y, max_z));
      // Ranges for border to be received.
      received_borders_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c, -halo_width_, 0    , 0),
             vec(c, -1          , max_y, max_z));
      received_borders_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c, max_x + 1           , 0    , 0),
             vec(c, max_x + halo_width_ , max_y, max_z));
      received_borders_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c, 0    , -halo_width_, 0),
             vec(c, max_x, -1          , max_z));
      received_borders_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c, 0    , max_y +1           , 0),
             vec(c, max_x, max_y + halo_width_, max_z));
      received_borders_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c, 0    , 0    , -halo_width_),
             vec(c, max_x, max_y, -1));
      received_borders_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c, 0    , 0    , max_z + 1),
             vec(c, max_x, max_y, max_z + halo_width_));
    }  // end of typedef block
    // Resize buffer for border exchange.
    borders_to_send_.resize(kDimensions*2);
    borders_to_send_(kBorderLeft).resize(number_of_components_to_exchange_,
                                         halo_width_, max_y + 1, max_z + 1);
    borders_to_send_(kBorderRight).resize(borders_to_send_(kBorderLeft)
                                          .shape());
    borders_to_send_(kBorderBottom).resize(number_of_components_to_exchange_,
                                           max_x + 1, halo_width_, max_z + 1);
    borders_to_send_(kBorderTop).resize(borders_to_send_(kBorderBottom)
                                        .shape());
    borders_to_send_(kBorderBack).resize(number_of_components_to_exchange_,
                                         max_x + 1, max_y + 1, halo_width_);
    borders_to_send_(kBorderFront).resize(borders_to_send_(kBorderBack)
                                          .shape());
    received_borders_.resize(kDimensions*2);
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      received_borders_(border).resize(borders_to_send_(border).shape());
      borders_to_send_(border) = 0;
      received_borders_(border) = 0;
      if (!borders_to_send_.isStorageContiguous() ||
          !received_borders_.isStorageContiguous()) {
        printf("Proc[%i]: Error! Halo exchange buffer is not contiguous!\n",
               process_rank_);
	throw std::invalid_argument("Error! Halo exchange buffer is not contiguous!\n");
      }  // end of if error
    }  // end of for border
    status_ = kSimulationStatusInitiated;
    return kDone;
  }  // end of BasicSimulationCore::Init()

  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Auto set boundary conditions for reduced dimenstions.
  ///
  /// @see BasicSimulationCore::Init() data_snapshot_ resizing.
  int SimulationInputConfig::AutoSetReducedBoundaryConditions() {
    if (grid_input_config_.get_total_grid_length(kAxisX) == 1) {
      boundary_condition_[kBorderRight] = kBoundaryConditionReduced;
      boundary_condition_[kBorderLeft] = kBoundaryConditionReduced;
    }  // end of if kAxisX dimension is reduced.
    if (grid_input_config_.get_total_grid_length(kAxisY) == 1) {
      boundary_condition_[kBorderTop] = kBoundaryConditionReduced;
      boundary_condition_[kBorderBottom] = kBoundaryConditionReduced;
    }  // end of if kAxisY dimension is reduced.
    if (grid_input_config_.get_total_grid_length(kAxisZ) == 1) {
      boundary_condition_[kBorderFront] = kBoundaryConditionReduced;
      boundary_condition_[kBorderBack] = kBoundaryConditionReduced;
    }  // end of if kAxisZ dimension is reduced.
    return kDone;
  }  // end of SimulationInputConfig::AutoSetReducedBoundaryConditions()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Parse config file keys.
  ///
  /// Read from config files keys (pair of name and value) and write them to
  /// config_file_map_ private member.
  int SimulationInputConfig::SetConfigFileMap() {
    FILE *config_file_pointer = fopen(config_file_name_.c_str(), "r");
    if (config_file_pointer == NULL) 
      throw std::invalid_argument("Error! Was not able to open config file!\n");
    std::string key_name, key_value;
    int reading_key_name = 1;
    config_file_map_.clear();
    while (!feof(config_file_pointer)) {
      int input_char = fgetc(config_file_pointer);
      if (input_char == '#') {  // Ignore comments till end of the line.
        while (input_char != '\n' && input_char != EOF)
          input_char = fgetc(config_file_pointer);
      }  // end of if comment needs to be ignored
      if (input_char == ' ') continue;  // Ignore spaces in config file.
      if (input_char == '=') {
        reading_key_name = 0;  // Stop reading key name.
        continue;
      }  // end of if start reading key value
      if (input_char == EOF || input_char == '\n') {
        reading_key_name = 1;
        if (key_name.size() != 0 && key_value.size() != 0) {
          if (config_file_map_.count(key_name)) 
	    throw std::invalid_argument
	      ("Error! Redefining key in config file is not allowed!\n");
          config_file_map_[key_name] = key_value;
        }  // end of if input is valid
        key_name = "";
        key_value = "";
        continue;
      }  // end of if end of line or end of file
      if (reading_key_name)
        key_name.push_back(input_char);
      else
        key_value.push_back(input_char);
    }  // end of while not end of config file
    fclose(config_file_pointer);
    return kDone;
  }  // end of SimulationInputConfig::SetConfigFileMap()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Check total PML width to be less than domain width.
  ///
  ///
  int SimulationInputConfig::CheckTotalPMLWidth() {
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      int total_pml_width = 0;
      if (boundary_condition_[axis] == kBoundaryConditionPML)
        total_pml_width += pml_width_;
      if (boundary_condition_[axis+kDimensions] == kBoundaryConditionPML)
        total_pml_width += pml_width_;
      if (total_pml_width >= grid_input_config_.
          get_total_grid_length(static_cast<Axis>(axis))) {
        status_ = kInputConfigErrorWidePml;
        int process_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
        if (process_rank == 0)
	  throw std::invalid_argument("Error! Too wide PML in axis = "
				      + std::to_string(axis));
      }  // end of if pml is too wide
    }  // end of for axis check pml width
    return kDone;
  }  // end of SimulationInputConfig::CheckTotalPMLWidth()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Set boundary conditions in all directions to be PML.
  ///
  ///
  int SimulationInputConfig::SetBoundaryConditionsAllPML() {
    for (int border = kBorderLeft; border < kDimensions*2; ++border)
      boundary_condition_[border] = kBoundaryConditionPML;
    return kDone;
  }  // end of SimulationInputConfig::SetBoundaryConditionsAllPML()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Preset config 1Dzero.
  ///
  /// At time step 231 field Ez is very close to zero.
  int SimulationInputConfig::PresetX1Dzero(int64_t total_time_steps,
                                          int64_t length_x,
                                          int64_t length_y,
                                          int64_t length_z) {
    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    if (process_rank == kOutput)
      printf("Preset X1Dzero\n");
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    SetBoundaryConditionsAllPML();
    pml_width_ = 1;
    pml_computational_ratio_ = 1.0;
    snapshot_interval_ = 1;
    halo_width_ = 1;
    time_depth_ = 2;
    algorithm_ = kAlgorithmSimpleX1D;
    // For most simple case of 1D  we will need Ez, Hy, epsilon, srcEz.
    number_of_grid_data_components_ = 4;
    number_of_components_to_exchange_ = 2;
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = kEz, kHy;
    total_time_steps_ = total_time_steps;
    return kDone;
  }  // end of SimulationInputConfig::PresetX1Dzero()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Preset config 2Dspeedup.
  ///
  /// Check speedup in cluster enviroment. Should be very nice for 16 nodes.
  /// @todo2(tig) Check preset
  int SimulationInputConfig::PresetTMz2Dspeedup(int64_t total_time_steps,
                                             int64_t length_x,
                                             int64_t length_y,
                                             int64_t length_z) {
    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    if (process_rank == kOutput)
      printf("Preset TMz2Dspeedup\n");
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    SetBoundaryConditionsAllPML();
    pml_width_ = 1;
    pml_computational_ratio_ = 1.0;
    snapshot_interval_ = 5;
    halo_width_ = 1;
    time_depth_ = 2;
    algorithm_ = kAlgorithmSimpleTMz2D;
    number_of_grid_data_components_ = 10;
    number_of_components_to_exchange_ = 3;
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = kEz, kHy, kHx;
    total_time_steps_ = total_time_steps;
    return kDone;
  }  // end of SimulationInputConfig::PresetTMz2Dspeedup()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Preset config 3Dsimple.
  ///
  /// Simply check 3D is runing.
  /// @todo2(tig) Check preset
  int SimulationInputConfig::Preset3Dsimple(int64_t total_time_steps,
                                            int64_t length_x,
                                            int64_t length_y,
                                            int64_t length_z) {
    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    if (process_rank == kOutput)
      printf("Preset 3Dsimple\n");
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    SetBoundaryConditionsAllPML();
    pml_width_ = 1;
    pml_computational_ratio_ = 1.0;
    snapshot_interval_ = 5;
    halo_width_ = 1;
    time_depth_ = 2;
    algorithm_ = kAlgorithmSimple3D;
    number_of_grid_data_components_ = 20;
    number_of_components_to_exchange_ = 6;
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = kEx, kEy, kEz, kHx, kHy, kHz;
    total_time_steps_ = total_time_steps;
    return kDone;
  }  // end of SimulationInputConfig::Preset3Dsimple()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Read top level config file
  ///
  /// @todo3 Currently values to read from config file are hard coded
  /// in ReadConfig(). Read them from real config file. Return some error
  /// for case if config file couldn be read.
  ///
  /// @todo1 Compare single core with meep. Current (rev40) state is
  /// not very good
  ///
  /// onza 128 proc timestep 100 size 6000x6000x1 time 11s.
  /// meep 32 processes 15 x 15 mkm res400
  /// on time step 235169 (time=293.961), 0.0897154 s/step
  /// Field time usage:
  /// connnecting chunks: 2.5294 s
  ///      time stepping: 20802.1 s
  ///      communicating: 439.32 s
  ///  outputting fields: 1.53428 s
  /// Fourier transforming: 0.128974 s
  ///    everything else: 60.4937 s
  ///
  /// Length of whole model
  /// 1 x 16 000 x 16 000 vertices x 8 components = 16 Gb on deb00
  /// 630 x 630 x 630 vertices x 8 components = 16 Gb on deb00
  int SimulationInputConfig::ReadConfig() {
    // Read config from file to private member.
    SetConfigFileMap();    
    int64_t length_x = 1, length_y = 1, length_z = 1;
    int64_t total_time_steps = 1;
    if (config_file_map_.count("length_x"))
      if (!convert_to_positive(config_file_map_["length_x"], length_x))
	throw std::invalid_argument
	  (std::string("Error! Non positive or non decimal input ") +
	   +"length_x = " + config_file_map_["length_x"] + "!\n" );
    if (config_file_map_.count("length_y"))
      if (!convert_to_positive(config_file_map_["length_y"], length_y))
    	throw std::invalid_argument
    	  (std::string("Error! Non positive or non decimal input ") +
    	   +"length_y = " + config_file_map_["length_y"] + "!\n" );
    if (config_file_map_.count("length_z"))
      if (!convert_to_positive(config_file_map_["length_z"], length_z))
	throw std::invalid_argument
	  (std::string("Error! Non positive or non decimal input ") +
	   +"length_z = " + config_file_map_["length_z"] + "!\n" );
    if (config_file_map_.count("total_time_steps"))
      if (!convert_to_positive(config_file_map_["total_time_steps"],
                               total_time_steps))
	throw std::invalid_argument
	  (std::string("Error! Non positive or non decimal input ")
	   + "total_time_steps" + config_file_map_["total_time_steps"] + "!\n");
    // convert_to_positive(config_file_map_["total_time_steps"]);
    // Select from some of hard coded tests.
    if (config_file_map_.count("test_case")) {
      if (config_file_map_["test_case"] == "X1Dzero") {
        PresetX1Dzero(total_time_steps, length_x, length_y, length_z);
      } else if (config_file_map_["test_case"] == "TMz2Dspeedup") {
	PresetTMz2Dspeedup(total_time_steps, length_x, length_y, length_z);
      } else if (config_file_map_["test_case"] == "3Dsimple") {
        Preset3Dsimple(total_time_steps, length_x, length_y, length_z);
      }
    } else {  // end of if simulating test case
      SetConfigFromFileMap();
    }  // end of else using castom settings.
    /// @todo1 For CPML implementation see Taflove 3d ed. p.307 section 7.9.2
    // Config auto check.
    AutoSetReducedBoundaryConditions();
    CheckTotalPMLWidth();
    status_ = kInputConfigAllDone;
    return kDone;
  };  // end of SimulationInputConfig::ReadConfig
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Set simulation parameters defined in config file.
  ///
  /// If presets are not used - parse parameters, validate their
  /// correctness, set values to SimulationInputConfig members.
  int SimulationInputConfig::SetConfigFromFileMap() {
    // todo2(tig) Do function realization.
    int process_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    if (process_rank == kOutput)
      printf("Using custom config\n");
    int64_t length_x = 200, length_y = 1, length_z = 1;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    // Setting boundary_condition_
    //            kBoundaryConditionPML or kBoundaryConditionPeriodical.
    SetBoundaryConditionsAllPML();
    pml_width_ = 1;
    pml_computational_ratio_ = 1.0;
    snapshot_interval_ = 5;
    halo_width_ = 1;
    time_depth_ = 2;
    algorithm_ = kAlgorithmSimpleX1D;
    number_of_grid_data_components_ = 4;
    int components_to_exchange[] = {kEz, kHy};
    total_time_steps_ = 240;  // check zero in 1D for Ez;
    number_of_components_to_exchange_ = sizeof(components_to_exchange)
                                        /sizeof(components_to_exchange[0]);
    components_to_exchange_.resize(number_of_components_to_exchange_);
    for (int i = 0; i < number_of_components_to_exchange_; ++i)
      components_to_exchange_(i) = components_to_exchange[i];
    //end of debug
    return kDone;
  }  // end of SimulationInputConfig::SetConfigFromFileMap()
}  // end of namespace onza
