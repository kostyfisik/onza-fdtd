///
/// @file   basic-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed May  2 13:42:32 2012
///
/// @brief  Simulation data and stepping algorithm
#include <mpi.h>
#include <blitz/array.h>
#include <iostream>
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
    local_time_step_++;
    if (total_time_steps_ < local_time_step_) return  kSimulationStatusFinished;
    return kSimulationStatusRunning;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Snapshot component from data_ to fiel.
  ///
  /// @warning  It is repeated MANY times! 
  void BasicSimulationCore::Snapshot() {
    if ( snapshot_interval_ * snapshot_frame_ < local_time_step_) {
      snapshot_frame_++;
      // if (total_time_steps_ - local_time_step_ > 20*snapshot_interval_) return;
      char filename[100];
      sprintf(filename,"Ez(x)[%i]-time%.7li-domain(x%.3dy%.3dz%.3d).onza",
              algorithm_, local_time_step_,
              my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ]);      
      FILE *snapshot;
      snapshot = fopen(filename, "w");
      int offset_x = my_coords_[kAxisX] * subdomain_size_[kAxisX];
      /// @todo2 make snapshot output  call be warper of algorithm
      /// dependent function call (call snapshot  with function pointer
      /// assigned at init time)
      switch (algorithm_) {
      case kAlgorithmSimpleX1D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot,"%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  data_(static_cast<int>(kEz), i, 0, 0));
        break;
      case kAlgorithmSimpleY1D:
        for (int i = 0; i < subdomain_size_[kAxisY]; ++i)
          fprintf(snapshot,"%li\t%g\n", i+subdomain_start_index_[kAxisY],
                  data_(static_cast<int>(kEz), 0, i, 0));
        break;
      case kAlgorithmSimpleZ1D:
        for (int i = 0; i < subdomain_size_[kAxisZ]; ++i)
          fprintf(snapshot,"%li\t%g\n", i+subdomain_start_index_[kAxisZ],
                  data_(static_cast<int>(kEz), 0, 0, i));
        break;
      case kAlgorithmSimpleTMz2D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot,"%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  data_(static_cast<int>(kEz), i, 140, 0));        
        break;
      case kAlgorithmSimple3D:
        for (int i = 0; i < subdomain_size_[kAxisX]; ++i)
          fprintf(snapshot,"%li\t%g\n", i+subdomain_start_index_[kAxisX],
                  // data_(static_cast<int>(kEz), i, 10, 10)); // 20x20
                  data_(static_cast<int>(kEz), i, 30, 30)); // 60x60
        break;
        // default:
        //   printf("Error! Should use some FDTD algorithm!\n");
        //   return kErrorWrongAlgorithm;
      }  // end of switch algorithm
      fclose(snapshot);
    }  // end of if it is time to output snapshot
  }  //end of BasicSimulationCore::Snapshot()
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
    if (process_rank_ == 0)
      switch (algorithm_) {
      case kAlgorithmSimpleX1D:
        data_(kSrcEz, 50, all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleY1D:
        data_(kSrcEz, all_x_, 50, all_z_) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleZ1D:
        data_(kSrcEz, all_x_, all_y_, 50) =
          exp(-pow2(local_time_step_ - 30.) / 100.);
        break;
      case kAlgorithmSimpleTMz2D:
        // data_(static_cast<int>(kSrcEz), 110, 140, all_z_) =
        //   exp(-pow2(local_time_step_ - 50.) / 200.);        
        data_(static_cast<int>(kSrcEz), 110, all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 50.) / 200.);        
        break;
      case kAlgorithmSimple3D:
        //point
        // data_(static_cast<int>(kSrcEz), 50, 30, 30) =
        //   exp(-pow2(local_time_step_ - 300.) / 4000.);
        //line
        // data_(static_cast<int>(kSrcEz), 50, 10, all_z_) =
        //   exp(-pow2(local_time_step_ - 100.) / 400.);
        // plane
        data_(static_cast<int>(kSrcEz), 50, all_y_, all_z_) =
          exp(-pow2(local_time_step_ - 100.) / 400.);
        break;
      // default:
      //   printf("Error! Should use some FDTD algorithm!\n");
      //   return kErrorWrongAlgorithm;
      }  // end of switch algorithm
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
    int c = 0;
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    // // Update H
    // // Hx
    // data_snapshot_        (t+1)(kHx  , x    , y    , z    ) =
    //   data_snapshot_      (t  )(kChxh, x    , y    , z    )
    //   * data_snapshot_    (t  )(kHx  , x    , y    , z    )
    //   +
    //   data_snapshot_      (t  )(kChxe, x    , y    , z    )
    //   * ((data_snapshot_  (t  )(kEy  , x    , y    , z + 1)
    //       - data_snapshot_(t  )(kEy  , x    , y    , z    ))
    //      -
    //      (data_snapshot_  (t  )(kEz  , x    , y + 1, z    )
    //       - data_snapshot_(t  )(kEz  , x    , y    , z    )));
    // // Hy
    // data_snapshot_        (t+1)(kHy  , x    , y    , z    ) =
    //   data_snapshot_      (t  )(kChyh, x    , y    , z    )
    //   * data_snapshot_    (t  )(kHy  , x    , y    , z    )
    //   +
    //   data_snapshot_      (t  )(kChye, x    , y    , z    )
    //   * ((data_snapshot_  (t  )(kEz  , x + 1, y    , z    )
    //       - data_snapshot_(t  )(kEz  , x    , y    , z    ))
    //      -
    //      (data_snapshot_  (t  )(kEx  , x    , y    , z + 1)
    //       - data_snapshot_(t  )(kEx  , x    , y    , z    )));
    // // Hz
    // data_snapshot_        (t+1)(kHz  , x    , y    , z    ) =
    //   data_snapshot_      (t  )(kChzh, x    , y    , z    )
    //   * data_snapshot_    (t  )(kHz  , x    , y    , z    )
    //   +
    //   data_snapshot_      (t  )(kChze, x    , y    , z    )
    //   * ((data_snapshot_  (t  )(kEx  , x    , y + 1, z    )
    //       - data_snapshot_(t  )(kEx  , x    , y    , z    ))
    //      -
    //      (data_snapshot_  (t  )(kEy  , x + 1, y    , z    )
    //       - data_snapshot_(t  )(kEy  , x    , y    , z    )));
    // // Update E
    // // Ex
    // c = 1;
    // data_snapshot_        (t+1)(kEx  , x    , y     + c, z     + c) =
    //   data_snapshot_      (t  )(kCexe, x    , y     + c, z     + c)
    //   * data_snapshot_    (t  )(kEx  , x    , y     + c, z     + c)
    //   +
    //   data_snapshot_      (t  )(kCexh, x    , y     + c, z     + c)
    //   * ((data_snapshot_  (t+1)(kHz  , x    , y     + c, z     + c)
    //       - data_snapshot_(t+1)(kHz  , x    , y - 1 + c, z     + c))
    //      -
    //      (data_snapshot_  (t+1)(kHy  , x    , y     + c, z     + c)
    //       - data_snapshot_(t+1)(kHy  , x    , y     + c, z - 1 + c)));
    // // Ey
    // data_snapshot_        (t+1)(kEy  , x     + c, y    , z     + c) =
    //   data_snapshot_      (t  )(kCeye, x     + c, y    , z     + c)
    //   * data_snapshot_    (t  )(kEy  , x     + c, y    , z     + c)
    //   +
    //   data_snapshot_      (t  )(kCeyh, x     + c, y    , z     + c)
    //   * ((data_snapshot_  (t+1)(kHx  , x     + c, y    , z     + c)
    //       - data_snapshot_(t+1)(kHx  , x     + c, y    , z - 1 + c))
    //      -
    //      (data_snapshot_  (t+1)(kHz  , x     + c, y    , z     + c)
    //       - data_snapshot_(t+1)(kHz  , x - 1 + c, y    , z     + c)));
    // // Ez
    // data_snapshot_        (t+1)(kEz   , x     + c, y     + c, z    ) =
    //   data_snapshot_      (t  )(kCeze , x     + c, y     + c, z    )
    //   * data_snapshot_    (t  )(kEz   , x     + c, y     + c, z    )
    //   +
    //   data_snapshot_      (t  )(kCezh , x     + c, y     + c, z    )
    //   * ((data_snapshot_  (t+1)(kHy   , x     + c, y     + c, z    )
    //       - data_snapshot_(t+1)(kHy   , x - 1 + c, y     + c, z    ))
    //      -
    //      (data_snapshot_  (t+1)(kHx   , x     + c, y     + c, z    )
    //       - data_snapshot_(t+1)(kHx   , x     + c, y - 1 + c, z    )))
    //   +
    //   data_snapshot_      (t  )(kSrcEz, x     + c, y     + c, z    );
    //debug Update H
    //debug Hx
    data_snapshot_        (t+1)(kHx  , x    , y    , z    ) =
      data_snapshot_      (t  )(kChxh, x    , y    , z    )
      * data_snapshot_    (t  )(kHx  , x    , y    , z    )
      +
      data_snapshot_      (t  )(kChxe, x    , y    , z    )
      * ((data_snapshot_  (t  )(kEy  , x    , y    , z + 1)
          - data_snapshot_(t  )(kEy  , x    , y    , z    ))
         -
         (data_snapshot_  (t  )(kEz  , x    , y + 1, z    )
          - data_snapshot_(t  )(kEz  , x    , y    , z    )));
    //debug Hy
    data_snapshot_        (t+1)(kHy  , x    , y    , z    ) =
      data_snapshot_      (t  )(kChyh, x    , y    , z    )
      * data_snapshot_    (t  )(kHy  , x    , y    , z    )
      +
      data_snapshot_      (t  )(kChye, x    , y    , z    )
      * ((data_snapshot_  (t  )(kEz  , x + 1, y    , z    )
          - data_snapshot_(t  )(kEz  , x    , y    , z    ))
         -
         (data_snapshot_  (t  )(kEx  , x    , y    , z + 1)
          - data_snapshot_(t  )(kEx  , x    , y    , z    )));
    //debug Hz
    data_snapshot_        (t+1)(kHz  , x    , y    , z    ) =
      data_snapshot_      (t  )(kChzh, x    , y    , z    )
      * data_snapshot_    (t  )(kHz  , x    , y    , z    )
      +
      data_snapshot_      (t  )(kChze, x    , y    , z    )
      * ((data_snapshot_  (t  )(kEx  , x    , y + 1, z    )
          - data_snapshot_(t  )(kEx  , x    , y    , z    ))
         -
         (data_snapshot_  (t  )(kEy  , x + 1, y    , z    )
          - data_snapshot_(t  )(kEy  , x    , y    , z    )));
    // Update E
    //debug Ex
    c = 1;
    data_snapshot_        (t+1)(kEx  , x    , y     + c, z     + c) =
      data_snapshot_      (t  )(kCexe, x    , y     + c, z     + c)
      * data_snapshot_    (t  )(kEx  , x    , y     + c, z     + c)
      +
      data_snapshot_      (t  )(kCexh, x    , y     + c, z     + c)
      * ((data_snapshot_  (t+1)(kHz  , x    , y     + c, z     + c)
          - data_snapshot_(t+1)(kHz  , x    , y - 1 + c, z     + c))
         -
         (data_snapshot_  (t+1)(kHy  , x    , y     + c, z     + c)
          - data_snapshot_(t+1)(kHy  , x    , y     + c, z - 1 + c)));
    //debug Ey
    data_snapshot_        (t+1)(kEy  , x     + c, y    , z     + c) =
      data_snapshot_      (t  )(kCeye, x     + c, y    , z     + c)
      * data_snapshot_    (t  )(kEy  , x     + c, y    , z     + c)
      +
      data_snapshot_      (t  )(kCeyh, x     + c, y    , z     + c)
      * ((data_snapshot_  (t+1)(kHx  , x     + c, y    , z     + c)
          - data_snapshot_(t+1)(kHx  , x     + c, y    , z - 1 + c))
         -
         (data_snapshot_  (t+1)(kHz  , x     + c, y    , z     + c)
          - data_snapshot_(t+1)(kHz  , x - 1 + c, y    , z     + c)));
    //debug Ez
    data_snapshot_        (t+1)(kEz   , x     + c, y     + c, z    ) =
      data_snapshot_      (t  )(kCeze , x     + c, y     + c, z    )
      * data_snapshot_    (t  )(kEz   , x     + c, y     + c, z    )
      +
      data_snapshot_      (t  )(kCezh , x     + c, y     + c, z    )
      * ((data_snapshot_  (t+1)(kHy   , x     + c, y     + c, z    )
          - data_snapshot_(t+1)(kHy   , x - 1 + c, y     + c, z    ))
         -
         (data_snapshot_  (t+1)(kHx   , x     + c, y     + c, z    )
          - data_snapshot_(t+1)(kHx   , x     + c, y - 1 + c, z    )))
      +
      data_snapshot_      (t  )(kSrcEz, x     + c, y     + c, z    );

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
    int c = 0;
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    // Update H
    // Hx
    data_snapshot_        (t+1)(kHx  , x    , y    , z    ) =
      data_snapshot_      (t  )(kChxh, x    , y    , z    )
      * data_snapshot_    (t  )(kHx  , x    , y    , z    )
      -
      data_snapshot_      (t  )(kChxe, x    , y    , z    )
      * (data_snapshot_   (t  )(kEz  , x    , y + 1, z    )
         - data_snapshot_ (t  )(kEz  , x    , y    , z    ));
    // Hy
    data_snapshot_        (t+1)(kHy  , x    , y    , z    ) =
      data_snapshot_      (t  )(kChyh, x    , y    , z    )
      * data_snapshot_    (t  )(kHy  , x    , y    , z    )
      +
      data_snapshot_      (t  )(kChye, x    , y    , z    )
      * (data_snapshot_   (t  )(kEz  , x + 1, y    , z    )
         - data_snapshot_ (t  )(kEz  , x    , y    , z    ));
    // Update E
    // Ez
    c=1;
    data_snapshot_        (t+1)(kEz   , x     + c, y     + c, z    ) =
      data_snapshot_      (t  )(kCeze , x     + c, y     + c, z    )
      * data_snapshot_    (t  )(kEz   , x     + c, y     + c, z    )
      +
      data_snapshot_      (t  )(kCezh , x     + c, y     + c, z    )
      * ((data_snapshot_  (t+1)(kHy   , x     + c, y     + c, z    )
          - data_snapshot_(t+1)(kHy   , x - 1 + c, y     + c, z    ))
         -
         (data_snapshot_  (t+1)(kHx   , x     + c, y     + c, z    )
          - data_snapshot_(t+1)(kHx   , x     + c, y - 1 + c, z    )))
      +
      data_snapshot_      (t  )(kSrcEz, x     + c, y     + c, z    );
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
    int c = 0;
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (t+1)(kHy, x     + c, y, z) =
      data_snapshot_  (t  )(kHy, x     + c, y, z)
      + data_snapshot_(t  )(kEz, x + 1 + c, y, z) * inv_imp0
      - data_snapshot_(t  )(kEz, x     + c, y, z) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (t+1)(kEz    , x     + c, y, z) =
        data_snapshot_  (t  )(kEz    , x     + c, y, z)
      + data_snapshot_  (t+1)(kHy    , x     + c, y, z)
        * data_snapshot_(t  )(kInvEps, x     + c, y, z)
        * imp0
      - data_snapshot_  (t+1)(kHy    , x - 1 + c, y, z)
        * data_snapshot_(t  )(kInvEps, x     + c, y, z)
        * imp0
      + data_snapshot_  (t  )(kSrcEz , x     + c, y, z);
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
    int c = 0;
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (t+1)(kHy, x, y     + c, z) =
      data_snapshot_  (t  )(kHy, x, y     + c, z)
      + data_snapshot_(t  )(kEz, x, y + 1 + c, z) * inv_imp0
      - data_snapshot_(t  )(kEz, x, y     + c, z) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (t+1)(kEz    , x, y     + c, z) =
        data_snapshot_  (t  )(kEz    , x, y     + c, z)
      + data_snapshot_  (t+1)(kHy    , x, y     + c, z)
        * data_snapshot_(t  )(kInvEps, x, y     + c, z)
        * imp0
      - data_snapshot_  (t+1)(kHy    , x, y - 1 + c, z)
        * data_snapshot_(t  )(kInvEps, x, y     + c, z)
        * imp0
      + data_snapshot_  (t  )(kSrcEz , x, y     + c, z);
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
    int c = 0;
    // Hy time = -1/2, x = +1/2, Ez time = 0, x = 0
    data_snapshot_    (t+1)(kHy, x, y, z     + c) =
      data_snapshot_  (t  )(kHy, x, y, z     + c)
      + data_snapshot_(t  )(kEz, x, y, z + 1 + c) * inv_imp0
      - data_snapshot_(t  )(kEz, x, y, z     + c) * inv_imp0;
    // Ez
    c = 1;
    data_snapshot_      (t+1)(kEz    , x, y, z     + c) =
        data_snapshot_  (t  )(kEz    , x, y, z     + c)
      + data_snapshot_  (t+1)(kHy    , x, y, z     + c)
        * data_snapshot_(t  )(kInvEps, x, y, z     + c)
        * imp0
      - data_snapshot_  (t+1)(kHy    , x, y, z - 1 + c)
        * data_snapshot_(t  )(kInvEps, x, y, z     + c)
        * imp0
      + data_snapshot_  (t  )(kSrcEz , x, y, z     + c);
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
    //debug
    // if (process_rank_ == kOutput && local_time_step_ >= 50 &&
    //     local_time_step_ < 56) {
    //   int z = 0;
    //   using namespace blitz;
    //   std::cout << "Cycle Do step["<<local_time_step_ <<"}: "
    //             << data_snapshot_(time_depth_-2)(Range(0,1),
    //                                              all_x_,
    //                                              0, 0) << std::endl;
    // }//end of if debug
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
        int opposite_border = (border + kDimensions) % (kDimensions*2);
        if (neighbours_ranks_[border] != MPI_PROC_NULL)
          data_(received_borders_range_(border, component))
                  (the_only_component_index, all_x_, all_y_, all_z_)
            = received_borders_(border)(component, all_x_, all_y_, all_z_);
      }  // end of for component
    }  // end of for border
    //debug
    // if (process_rank_ == kOutput) {
    //   printf("Do border step..\n");
    //   std::cout<<"all_x :"<<all_data_x_ <<" all_y_:"<<all_data_y_<<std::endl;
    // }
    if (neighbours_ranks_[kBorderLeft] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(data_border_range_[kBorderLeft], inner_y_, inner_z_);
    //debug
    // if (process_rank_ == kOutput) printf("Do border step2..\n");
    if (neighbours_ranks_[kBorderRight] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(data_border_range_[kBorderRight], inner_y_, inner_z_);
    if (neighbours_ranks_[kBorderBottom] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(inner_x_, data_border_range_[kBorderBottom], inner_z_);
    if (neighbours_ranks_[kBorderTop] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(inner_x_, data_border_range_[kBorderTop], inner_z_);
    if (neighbours_ranks_[kBorderBack] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(inner_x_, inner_y_, data_border_range_[kBorderBack]);
    if (neighbours_ranks_[kBorderFront] != MPI_PROC_NULL)
      (this->*RunAlgorithm)(inner_x_, inner_y_, data_border_range_[kBorderFront]);
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
      data_(kInvEps, all_x_, all_y_, all_z_) = 1.0/epsilon;// / (process_rank_ + 1);
      break;
    case kAlgorithmSimpleTMz2D:
      //debug
      data_(kInvEps, all_x_, all_y_, all_z_) = 1.0/epsilon;// / (process_rank_ + 1);
      data_(kCeze, all_x_, all_y_, all_z_) = 1.0;
      data_(kCezh, all_x_, all_y_, all_z_) = courant_2D * imp0;
      data_(kChxh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChxe, all_x_, all_y_, all_z_) = courant_2D / imp0;
      data_(kChyh, all_x_, all_y_, all_z_) = 1.0;
      data_(kChye, all_x_, all_y_, all_z_) = courant_2D / imp0;
      break;
    case kAlgorithmSimple3D:
      //debug
      data_(kInvEps, all_x_, all_y_, all_z_) = 1.0/epsilon;// / (process_rank_ + 1);
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
      printf("Error! Should use some FDTD algorithm!\n");
      return kErrorWrongAlgorithm;
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
                                 int64_t subdomain_finish_index[])
  {
    if (simulation_input_config_.status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;    
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
      printf("Error! Should use some FDTD algorithm!\n");
      return kErrorWrongAlgorithm;
    }  // end of switch algorithm
    halo_width_ = simulation_input_config_.halo_width();
    if (halo_width_ < 1) {
      printf("Error! FDTD algorithm`s should have some subdomain halo!\n");
      return kErrorWrongHaloWidth;
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
      printf("Error! FDTD algorithm`s time depth should be >= 2!\n");
      return kErrorWrongTimeDepth;
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
      // reading config boundary conditions kBoudaryConditionReduced
      // should be set in this dimension.
      data_snapshot_(time).resize(blitz::shape(number_of_grid_data_components_,
               max_x + 1 + halo_width_x * 2,   max_y + 1 + halo_width_y * 2,
               max_z + 1 + halo_width_z * 2));
      data_snapshot_(time).reindexSelf(blitz::shape(0, -halo_width_x,
                                 -halo_width_y, -halo_width_z));
    }  // end of for time of snapshot
    data_.resize(blitz::shape(number_of_grid_data_components_,
                              max_x + 1 + halo_width_x * 2,   max_y + 1 + halo_width_y * 2,
                              max_z + 1 + halo_width_z * 2));
    data_.reindexSelf(blitz::shape(0, -halo_width_x,
                                   -halo_width_y, -halo_width_z));
    // Ranges for all data spacial coords.
    all_x_ = blitz::Range::all();
    all_y_ = blitz::Range::all();
    all_z_ = blitz::Range::all();
    //debug
    // all_data_x_ = blitz::Range(-halo_width_x, max_x + halo_width_x);
    // all_data_y_ = blitz::Range(-halo_width_y, max_y + halo_width_y);
    // all_data_z_ = blitz::Range(-halo_width_z, max_z + halo_width_z);
    // Special reference to current time step data.
    data_.reference(data_snapshot_(time_depth_-2));
    // Specify ranges in each dimension for FDTD time stepping. For
    // ease of defining FDTD equations range can be shifted in forward
    // or backwrd direction for any decimal value less or equal
    // halo_width_. Note that aligment shift for each FDTD equation
    // should be applied to provide all ranges to be valid.
    inner_x_ = blitz::Range(0, max_x);
    inner_y_ = blitz::Range(0, max_y);
    inner_z_ = blitz::Range(0, max_z);
    // Specify data range for each border. Due to choice of
    // aligment shift in FDTD equation low and hi index ranges have
    // different number of elements.
    data_border_range_[kBorderLeft] = blitz::Range(-halo_width_x, halo_width_x);
    data_border_range_[kBorderBottom] = blitz::Range(-halo_width_y, halo_width_y);
    data_border_range_[kBorderBack] = blitz::Range(-halo_width_z, halo_width_z);
    data_border_range_[kBorderRight] = blitz::Range(max_x - halo_width_x, max_x);
    data_border_range_[kBorderTop] = blitz::Range(max_y - halo_width_y, max_y);
    data_border_range_[kBorderFront] = blitz::Range(max_z - halo_width_z, max_z);
    // //debug
    // if (process_rank_ == 0) {
    //   for (int border = 0; border < kDimensions*2; ++border)        
    //     std::cout<<"border["<<border<<"]: "<<data_border_range_[border]<<std::endl;
    //   std::cout<<"ds(t) shape: "<< data_snapshot_(0).shape()<<std::endl
    //            <<"lb : "<<data_snapshot_(0).lbound()<<std::endl
    //            <<"ub : "<<data_snapshot_(0).ubound()
    //            <<std::endl<<std::flush;
    // }  //end of debug
    // Specify borders ranges for halo exchange.
    // Used for quick selection in BasicSimulationCore::PrepareBordersToSend().
    borders_to_send_range_.resize(kDimensions*2, number_of_components_to_exchange_);
    received_borders_range_.resize(kDimensions*2, number_of_components_to_exchange_);
    for (int component = 0; component < number_of_components_to_exchange_; ++component) {
      // Try to shorten notation
      typedef blitz::RectDomain<1+kDimensions> rd;
      typedef blitz::TinyVector<int,4> vec;
      int c = components_to_exchange_(component);
      // Ranges for border to be send.
      borders_to_send_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c, 0            , 0    , 0    ),
             vec(c, halo_width_-1, max_y, max_z));
      borders_to_send_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c, max_x - (halo_width_-1), 0    , 0    ),
             vec(c, max_x                  , max_y, max_z));
      borders_to_send_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c, 0    , 0            , 0    ),
             vec(c, max_x, halo_width_-1, max_z));
      borders_to_send_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c, 0    , max_y - (halo_width_-1), 0    ),
             vec(c, max_x, max_y                  , max_z));
      borders_to_send_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c, 0    , 0    , 0            ),
             vec(c, max_x, max_y, halo_width_-1));
      borders_to_send_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c, 0    , 0    , max_z - (halo_width_-1)),
             vec(c, max_x, max_y, max_z                  ));
      // Ranges for border to be received.
      received_borders_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c, -halo_width_, 0    , 0    ),
             vec(c, -1          , max_y, max_z));
      received_borders_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c, max_x + 1           , 0    , 0    ),
             vec(c, max_x + halo_width_ , max_y, max_z));
      received_borders_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c, 0    , -halo_width_, 0    ),
             vec(c, max_x, -1          , max_z));
      received_borders_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c, 0    , max_y +1           , 0    ),
             vec(c, max_x, max_y + halo_width_, max_z));
      received_borders_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c, 0    , 0    , -halo_width_),
             vec(c, max_x, max_y, -1          ));
      received_borders_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c, 0    , 0    , max_z + 1          ),
             vec(c, max_x, max_y, max_z + halo_width_));
    }  // end of typedef block
    // Resize buffer for border exchange.
    borders_to_send_.resize(kDimensions*2);
    borders_to_send_(kBorderLeft).resize(number_of_components_to_exchange_,
                                         halo_width_, max_y + 1, max_z + 1);
    borders_to_send_(kBorderRight).resize(borders_to_send_(kBorderLeft).shape());
    borders_to_send_(kBorderBottom).resize(number_of_components_to_exchange_,
                                           max_x + 1, halo_width_, max_z + 1);
    borders_to_send_(kBorderTop).resize(borders_to_send_(kBorderBottom).shape());
    borders_to_send_(kBorderBack).resize(number_of_components_to_exchange_,
                                         max_x + 1, max_y + 1, halo_width_);
    borders_to_send_(kBorderFront).resize(borders_to_send_(kBorderBack).shape());
    received_borders_.resize(kDimensions*2);
    for (int border = kBorderLeft; border < kDimensions*2; ++border) {
      received_borders_(border).resize(borders_to_send_(border).shape());
      borders_to_send_(border) = 0;
      received_borders_(border) = 0;
      if (!borders_to_send_.isStorageContiguous() ||
          !received_borders_.isStorageContiguous()) {
        printf("Proc[%i]: Error! Halo exchange buffer is not contiguous!\n", process_rank_);
        return kErrorExchangeBufferIsNotContiguous;
      }  // end of if error
    }  // end of for border
    status_ = kSimulationStatusInitiated;
    return kDone;
  }  // end of BasicSimulationCore::Init()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Read top level config file
  ///
  /// @todo3 Currently values to read from config file are hard coded
  /// in ReadConfig(). Read them from real config file. Return some error
  /// for case if config file couldn be read.
  int SimulationInputConfig::ReadConfig() {
    // ********************************************************************** //
    // Length of whole model
    // 1 x 16 000 x 16 000 vertices x 8 components = 16 Gb on deb00
    // 630 x 630 x 630 vertices x 8 components = 16 Gb on deb00
    // int64_t length_x = 6000, length_y = 6000, length_z = 1;  
    // int64_t length_x = 130, length_y = 130, length_z = 130;  
    // int64_t length_x = 113, length_y = 120, length_z = 179;  // !!
    // int64_t length_x = 101, length_y = 307, length_z = 908; // !!
    // int64_t length_x = 813, length_y = 1, length_z = 79; // !!
    // int64_t length_x = 5, length_y = 9, length_z = 2; // Best to go with MPIsize = 3 
    // int64_t length_x = 800, length_y = 1, length_z = 1;
    // int64_t length_x = 200, length_y = 180, length_z = 180;  // check 3D algorithm point source
    int64_t length_x = 200, length_y = 20, length_z = 20;  // check 3D algorithm !!!?small?!!!
    // int64_t length_x = 501, length_y = 281, length_z = 1; // check 2D algorithm TMz
    // int64_t length_x = 200, length_y = 1, length_z = 1; // check 1D algorithm X
    // int64_t length_x = 1, length_y = 200, length_z = 1; // check 1D algorithm Y
    // int64_t length_x = 1, length_y = 1, length_z = 200; // check 1D algorithm Z
    // int64_t length_x = 4, length_y = 1, length_z = 1;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    // ********************************************************************** //
    // Setting boundary_condition_.
    boundary_condition_[kBorderRight] = kBoudaryConditionPML;
    boundary_condition_[kBorderLeft] = kBoudaryConditionPML;
    // boundary_condition_[kBorderRight] = kBoudaryConditionPeriodical;
    // boundary_condition_[kBorderLeft] = kBoudaryConditionPeriodical;
    boundary_condition_[kBorderTop] = kBoudaryConditionPML;
    boundary_condition_[kBorderBottom] = kBoudaryConditionPML;
    // boundary_condition_[kBorderTop] = kBoudaryConditionPeriodical;
    // boundary_condition_[kBorderBottom] = kBoudaryConditionPeriodical;
    boundary_condition_[kBorderFront] = kBoudaryConditionPML;
    boundary_condition_[kBorderBack] = kBoudaryConditionPML;
    // boundary_condition_[kBorderFront] = kBoudaryConditionPeriodical;
    // boundary_condition_[kBorderBack] = kBoudaryConditionPeriodical;
    // ********************************************************************** //
    // Auto set periodical boundary condition for reduced dimenstions.
    // See also BasicSimulationCore::Init() data_snapshot_ resizing.
    if (length_x == 1) {
      boundary_condition_[kBorderRight] = kBoudaryConditionReduced;
      boundary_condition_[kBorderLeft] = kBoudaryConditionReduced;
    }  // end of if kAxisX dimension is reduced.
    if (length_y == 1) {
      boundary_condition_[kBorderTop] = kBoudaryConditionReduced;
      boundary_condition_[kBorderBottom] = kBoudaryConditionReduced;
    }  // end of if kAxisY dimension is reduced.
    if (length_z == 1) {
      boundary_condition_[kBorderFront] = kBoudaryConditionReduced;
      boundary_condition_[kBorderBack] = kBoudaryConditionReduced;
    }  // end of if kAxisZ dimension is reduced.
    /// For CPML implementation see Taflove 3d ed. p.307 section 7.9.2
    pml_width_ = 1;
    pml_computational_ratio_ = 1.0;
    // pml_width_ = 7;
    // pml_computational_ratio_ = 1.27;
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      int total_pml_width = 0;
      if (boundary_condition_[axis] == kBoudaryConditionPML)
        total_pml_width += pml_width_;
      if (boundary_condition_[axis+kDimensions] == kBoudaryConditionPML)
        total_pml_width += pml_width_;
      if (total_pml_width >= grid_input_config_.
          get_total_grid_length(static_cast<Axis>(axis))) {
        status_ = kInputConfigErrorWidePml;
        int process_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
        if (process_rank == 0)
          printf("Error! Too wide PML in axis = %i!\n", axis);
        return kErrorTooWidePml;
      }  // end of if pml is too wide
    }  // end of for axis check pml width
    // onza 128 proc timestep 100 size 6000x6000x1 time 33s.
    // meep 32 processes 15 x 15 mkm res400 
    // on time step 235169 (time=293.961), 0.0897154 s/step
    // Field time usage:
    // connnecting chunks: 2.5294 s
    //      time stepping: 20802.1 s
    //      communicating: 439.32 s
    //  outputting fields: 1.53428 s
    // Fourier transforming: 0.128974 s
    //    everything else: 60.4937 s
    // For most simple case of 1D  we will need Ez, Hy, epsilon, srcEz.
    snapshot_interval_ = 5.0;
    halo_width_ = 1;
    time_depth_ = 2;
    // ********************************************************************** //
    // 1D section
    // ********************************************************************** //
    // algorithm_ = kAlgorithmSimpleZ1D;
    // number_of_grid_data_components_ = 4;
    // int components_to_exchange[] = {kEz, kHy};
    // ********************************************************************** //
    // 2D section
    // ********************************************************************** //
    // algorithm_ = kAlgorithmSimpleTMz2D;
    // number_of_grid_data_components_ = 11;
    // int components_to_exchange[] = {kEz, kHy, kHx};
    // ********************************************************************** //
    // 3D section
    // ********************************************************************** //
    algorithm_ = kAlgorithmSimple3D;
    number_of_grid_data_components_ = 20;
    int components_to_exchange[] = {kEx, kEy, kEz, kHx, kHy, kHz};
    // ********************************************************************** //
    number_of_components_to_exchange_ = sizeof(components_to_exchange) / sizeof(int);
    components_to_exchange_.resize(number_of_components_to_exchange_);
    for (int i = 0; i < number_of_components_to_exchange_; ++i)
      components_to_exchange_(i) = components_to_exchange[i];
    // total_time_steps_ = 450000;
    total_time_steps_ = 640;
    //total_time_steps_ = 240; //check zero in 1D for Ez;
    // total_time_steps_ = 100;
    // total_time_steps_ = 60;
    // total_time_steps_ = 4;
    status_ = kInputConfigAllDone;
    return kDone;
  };  // end of SimulationInputConfig::ReadConfig
}  // end of namespace onza
