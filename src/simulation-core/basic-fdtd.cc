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
              = data_(borders_range_(border, component))(the_only_component_index,
                                                         all_x_, all_y_, all_z_);
          else
            received_borders_(opposite_border)(component, all_x_, all_y_, all_z_)
              = data_(borders_range_(border, component))(the_only_component_index,
                                                         all_x_, all_y_, all_z_);
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
      char filename[100];
      sprintf(filename,"Ez(x)-time%.7li-domain(x%.3dy%.3dz%.3d).onza",
              local_time_step_,
              my_coords_[kAxisX], my_coords_[kAxisY], my_coords_[kAxisZ]);      
      FILE *snapshot;
      snapshot = fopen(filename, "w");
      int data_x_dimension = 1;
      int offset_x = my_coords_[kAxisX] * data_.extent(data_x_dimension);
      for (int i = 0; i < data_.extent(data_x_dimension); ++i)
        fprintf(snapshot,"%i\t%g\n", offset_x + i,
                data_(static_cast<int>(kEz), i, 0, 0));
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
    //debug 1D
    if (process_rank_ == 0)
      data_(kSrcEz, subdomain_size_[kAxisX]/2, all_y_, all_z_) =
        exp(-pow2(local_time_step_ - 30.) / 100.);
  }  // end of BasicSimulationCore::PrepareSource()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief FDTD algorithm for 1D case (length in kAxisX dimenstion).
  ///
  /// @warning  It is repeated MANY times!
  /// @todo1 Check not doing data copying during algorithm call.
  void BasicSimulationCore::FDTD_1D_axis_x(
      blitz::Array<double, 1+kDimensions> data,
      blitz::Range x,
      blitz::Range y,
      blitz::Range z) {
    // 1D test FDTD algorithm
    /// @todo1 Implement variable sized array for constants.
    // impedance of free space.    
    double imp0 = 377.0;
    double inv_imp0 = 1.0/imp0;
    data(kHy    , x    , y, z) =
        data(kHy, x    , y, z)
      + data(kEz, x + 1, y, z) * inv_imp0
      - data(kEz, x    , y, z) * inv_imp0;
    //  ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0 / epsR[mm];
    data(kEz          , x  , y, z) =
        data(kEz      , x  , y, z)
      + data(kHy      , x  , y, z)
        * data(kInvEps, x  , y, z)
        * imp0
      - data(kHy      , x-1, y, z)
        * data(kInvEps, x  , y, z)
        * imp0
      + data(kSrcEz   , x  , y, z);
  }  // end of BasicSimulationCore::FDTD_1D_axis_x()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for border part of grid data.
  ///
  /// @warning  It is repeated MANY times! 
  void BasicSimulationCore::DoBorderStep() {
    //debug Hardly selecting FDTD algorithm.
    // FDTD_1D_axis_x(data_, inner_x_, inner_y_, inner_z_);
  }  // end of BasicSimulationCore::DoStep()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for internal part of grid data.
  ///
  /// @warning  It is repeated MANY times! 
  void BasicSimulationCore::DoStep() {
    //debug Hardly selecting FDTD algorithm.
    FDTD_1D_axis_x(data_, inner_x_, inner_y_, inner_z_);
  }  // end of BasicSimulationCore::DoStep()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize grid data for all components
  int BasicSimulationCore::SetGridData() {
    // set all fields (Ex, Ey, Ez, Hx, Hy, Hz) to be zero. And everything else.
    data_ = 0;
    //debug
    data_(kInvEps, all_x_, all_y_, all_z_) = 1;// / (process_rank_ + 1);
    //debug
    // if (process_rank_ == 0) {
    //   int z = 0;
    //   std::cout << data_(kEz, all_x_, all_y_, 0) << std::endl;
    // }
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
                                 int my_coords[])
  {
    for (int border = kBorderLeft; border < kDimensions*2; ++border)
      neighbours_ranks_[border] = neighbours_ranks[border];
    status_ = kSimulationStatusUndefined;
    if (simulation_input_config_.status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;    
    // Init member in alphabetical order    
    all_x_ = blitz::Range::all();
    all_y_ = blitz::Range::all();
    all_z_ = blitz::Range::all();
    halo_width_ = simulation_input_config_.halo_width();
    local_time_step_ = 0;
    for (int axis = kAxisX; axis < kDimensions; ++axis) {
      my_coords_[axis] = my_coords[axis];
    }
    number_of_components_to_exchange_ = simulation_input_config_.
      number_of_components_to_exchange();
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = simulation_input_config_.components_to_exchange();
    number_of_grid_data_components_ = simulation_input_config_.
      number_of_grid_data_components();
    process_rank_ = process_rank;
    snapshot_interval_ = simulation_input_config_.snapshot_interval();
    snapshot_frame_ = 0;
      for (int axis = kAxisX; axis < kDimensions; ++axis)
      subdomain_size_[axis] = subdomain_size[axis];
    total_time_steps_ = simulation_input_config_.total_time_steps();
    // Resize grid data to fit subdomain of current MPI process.
    int max_x = subdomain_size_[kAxisX]-1;
    int max_y = subdomain_size_[kAxisY]-1;
    int max_z = subdomain_size_[kAxisZ]-1;
    data_.resize(number_of_grid_data_components_, max_x + 1, max_y + 1, max_z + 1);
    // Specify inner range in each dimension.
    if (max_x - halo_width_ > 0)
      inner_x_ = blitz::Range(halo_width_, max_x - halo_width_);
    else
      inner_x_ = blitz::Range(0          , max_x);
    if (max_y - halo_width_ > 0)
      inner_y_ = blitz::Range(halo_width_, max_y - halo_width_);
    else
      inner_y_ = blitz::Range(0          , max_y);
    if (max_z - halo_width_ > 0)
      inner_z_ = blitz::Range(halo_width_, max_z - halo_width_);
    else
      inner_z_ = blitz::Range(0          , max_z);
    // Define ranges for borders selection.
    borders_range_.resize(kDimensions*2, number_of_components_to_exchange_);
    for (int component = 0; component < number_of_components_to_exchange_; ++component) {
      // Describe indexes for each border.
      typedef blitz::RectDomain<1+kDimensions> rd;
      typedef blitz::TinyVector<int,4> vec;
      int c = components_to_exchange_(component);
      borders_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c, 0            , 0    , 0    ),
             vec(c, halo_width_-1, max_y, max_z));
      borders_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c, max_x - (halo_width_-1), 0    , 0    ),
             vec(c, max_x                  , max_y, max_z));
      borders_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c, 0    , 0            , 0    ),
             vec(c, max_x, halo_width_-1, max_z));
      borders_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c, 0    , max_y - (halo_width_-1), 0    ),
             vec(c, max_x, max_y                  , max_z));
      borders_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c, 0    , 0    , 0            ),
             vec(c, max_x, max_y, halo_width_-1));
      borders_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c, 0    , 0    , max_z - (halo_width_-1)),
             vec(c, max_x, max_y, max_z                  ));
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
    // int64_t length_x = 1000, length_y = 1000, length_z = 1;  
    // int64_t length_x = 113, length_y = 120, length_z = 179;  // !!
    // int64_t length_x = 101, length_y = 307, length_z = 908; // !!
    // int64_t length_x = 813, length_y = 1, length_z = 79; // !!
    // int64_t length_x = 5, length_y = 9, length_z = 2; // Best to go with MPIsize = 3 
    int64_t length_x = 200, length_y = 1, length_z = 1;
    // int64_t length_x = 4, length_y = 1, length_z = 1;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    // ********************************************************************** //
    // Setting boundary_condition_.
    // boundary_condition_[kBorderRight] = kBoudaryConditionPML;
    // boundary_condition_[kBorderLeft] = kBoudaryConditionPML;
    boundary_condition_[kBorderRight] = kBoudaryConditionPeriodical;
    boundary_condition_[kBorderLeft] = kBoudaryConditionPeriodical;
    boundary_condition_[kBorderTop] = kBoudaryConditionPML;
    boundary_condition_[kBorderBottom] = kBoudaryConditionPML;
    // boundary_condition_[kBorderTop] = kBoudaryConditionPeriodical;
    // boundary_condition_[kBorderBottom] = kBoudaryConditionPeriodical;
    // boundary_condition_[kBorderFront] = kBoudaryConditionPML;
    // boundary_condition_[kBorderBack] = kBoudaryConditionPML;
    boundary_condition_[kBorderFront] = kBoudaryConditionPeriodical;
    boundary_condition_[kBorderBack] = kBoudaryConditionPeriodical;
    // ********************************************************************** //
    // Auto set periodical boundary condition for reduced dimenstions./
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
    pml_width_ = 2;
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
    number_of_grid_data_components_ = 4;
    int components_to_exchange[] = {kEz, kHy};
    // int components_to_exchange[] = {kEx, kEy};
    number_of_components_to_exchange_ = sizeof(components_to_exchange) / sizeof(int);
    components_to_exchange_.resize(number_of_components_to_exchange_);
    for (int i = 0; i < number_of_components_to_exchange_; ++i)
      components_to_exchange_(i) = components_to_exchange[i];
    total_time_steps_ = 450;
    status_ = kInputConfigAllDone;
    return kDone;
  };  // end of SimulationInputConfig::ReadConfig
}  // end of namespace onza
