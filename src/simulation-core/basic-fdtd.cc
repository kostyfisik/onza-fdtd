///
/// @file   basic-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed May  2 13:42:32 2012
///
/// @brief  Simulation data and stepping algorithm
#include <mpi.h>
#include <blitz/array.h>
#include <cstdio>
#include "./basic-fdtd.h"
#include "../common.h"
namespace onza {
  using blitz::Range;
  /// @brief Initialize process variables
  ///
  /// @return 0 or error code.
  ///
  int  BasicSimulationCore::Init() {
    if (simulation_input_config_.ReadConfig() != kDone)
      return kErrorTooWidePml;
    if (simulation_input_config_.status() != kInputConfigAllDone)
      return kErrorUsingInputConfigTooEarly;
    /// @todo2 use length values from subdomain properties
    int64_t length_x = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisX);
    int64_t length_y = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisY);
    int64_t length_z = simulation_input_config_.grid_input_config_.
        get_total_grid_length(kAxisZ);
    // test_array_.resize(length_x, length_y, length_z);
    // test_array_(0, 0, Range::all()) = 0;
    return kDone;
  }  // end of BasicSimulationCore::Init
  /// @brief Read top level config file
  ///
  /// @todo3 Currently values to read from config file are hard coded
  /// in ReadConfig(). Read them from real config file. Return some error
  /// for case if config file couldn be read.
  int SimulationInputConfig::ReadConfig() {
    // Length of whole model
    int64_t length_x = 113, length_y = 120, length_z = 79;  // !!
    // int64_t length_x = 101, length_y = 307, length_z = 908; // !!
    // int64_t length_x = 813, length_y = 1, length_z = 79; // !!
    // int64_t length_x = 100, length_y = 2, length_z = 2;
    // int64_t length_x = 360, length_y = 1, length_z = 1;
    // int64_t length_x = 4, length_y = 1, length_z = 1;
    grid_input_config_.set_total_grid_length(length_x, length_y, length_z);
    /// For CPML implementation see Taflove 3d ed. p.307 section 7.9.2
    pml_width_ = 7;
    pml_computational_ratio_ = 1.27;
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
