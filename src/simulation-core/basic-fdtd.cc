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
  using blitz::Range;
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for internal part of grid data.
  int BasicSimulationCore::DoStep() {
    local_time_step_++;
    if (total_time_steps_ - local_time_step_ < 1) return  kSimulationStatusFinished;
    return kSimulationStatusRunning;
  }  // end of BasicSimulationCore::DoStep()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize grid data for all components
  int BasicSimulationCore::SetGridData() {
    Range all_x = Range::all(), all_y = Range::all(), all_z = Range::all();
    int first = 0;
    // set all fields (Ex, Ey, Ez, Hx, Hy, Hz) to be zero.
    for (int component = kEx; component <= kHz; ++component)
      data_(component, all_x, all_y, all_z) = 0;
    blitz::firstIndex component_position_x;
    blitz::secondIndex component_position_y;
    blitz::thirdIndex component_position_z;
    blitz::Array<double, kDimensions> ex = data_(kEx, all_x, all_y, all_z); // view
    /// @todo1 Remove debug assigment for field ex.
    ex = (component_position_y + component_position_x * 10
          + ex.columns() * process_rank_)
      * (component_position_z*2 - 1) * (-1);
    if (process_rank_ == 0)
      std::cout << data_(kEx, all_x, all_y, 0);
    // if (ex.isStorageContiguous())
    //   printf("ex for proc %i is contiguous!\n", process_rank_);
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  int  BasicSimulationCore::Init(const int64_t subdomain_size[], int process_rank) {
    status_ = kSimulationStatusUndefined;
    if (simulation_input_config_.status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;
    // Init member in alphabetical order
    halo_width_ = simulation_input_config_.halo_width();
    local_time_step_ = 0;
    number_of_grid_data_components_ = simulation_input_config_.
      number_of_grid_data_components();
    process_rank_ = process_rank;
    for (int axis = kAxisX; axis < kDimensions; ++axis)
      subdomain_size_[axis] = subdomain_size[axis];
    total_time_steps_ = simulation_input_config_.total_time_steps();
    // Resize grid data to fit subdomain of current MPI process.
    data_.resize(number_of_grid_data_components_,
                 subdomain_size_[kAxisX],
                 subdomain_size_[kAxisY],
                 subdomain_size_[kAxisZ]);
    status_ = kSimulationStatusInitiated;
    return kDone;
  }  // end of BasicSimulationCore::Init
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Read top level config file
  ///
  /// @todo3 Currently values to read from config file are hard coded
  /// in ReadConfig(). Read them from real config file. Return some error
  /// for case if config file couldn be read.
  int SimulationInputConfig::ReadConfig() {
    total_time_steps_ = 3;
    // For most simple case we will need Ex, Ey, Ez, epsilon, Hx, Hy,
    // Hz, mu.
    number_of_grid_data_components_ = 8;
    halo_width_ = 1;
    /// For CPML implementation see Taflove 3d ed. p.307 section 7.9.2
    // pml_width_ = 7;
    // pml_computational_ratio_ = 1.27;
    pml_width_ = 2;
    pml_computational_ratio_ = 1.02;
    // ********************************************************************** //
    // Length of whole model
    // 1 x 16 000 x 16 000 vertices x 8 components = 16 Gb on deb00
    // 630 x 630 x 630 vertices x 8 components = 16 Gb on deb00
    // int64_t length_x = 1, length_y = 642, length_z = 380;  
    // int64_t length_x = 113, length_y = 120, length_z = 179;  // !!
    // int64_t length_x = 101, length_y = 307, length_z = 908; // !!
    // int64_t length_x = 813, length_y = 1, length_z = 79; // !!
    int64_t length_x = 7, length_y = 9, length_z = 2;
    // int64_t length_x = 360, length_y = 1, length_z = 1;
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
      boundary_condition_[kBorderRight] = kBoudaryConditionPeriodical;
      boundary_condition_[kBorderLeft] = kBoudaryConditionPeriodical;
    }  // end of if kAxisX dimension is reduced.
    if (length_y == 1) {
      boundary_condition_[kBorderTop] = kBoudaryConditionPeriodical;
      boundary_condition_[kBorderBottom] = kBoudaryConditionPeriodical;
    }  // end of if kAxisY dimension is reduced.
    if (length_z == 1) {
      boundary_condition_[kBorderFront] = kBoudaryConditionPeriodical;
      boundary_condition_[kBorderBack] = kBoudaryConditionPeriodical;
    }  // end of if kAxisZ dimension is reduced.
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
      }
    }
    status_ = kInputConfigAllDone;
    return kDone;
  };  // end of SimulationInputConfig::ReadConfig
}  // end of namespace onza
