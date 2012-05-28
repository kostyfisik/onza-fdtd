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
        // Get address of data_ slice for current border and component
        // In fist (component) dimension it will have the only element
        // (current component)
        blitz::Array<double, 1+kDimensions> reference_to_component_border_data
          = data_(borders_range_(border, component));
        // Get address of borders_to_send_ element for current border
        blitz::Array<double, 1+kDimensions> reference_to_component_to_send
          = borders_to_send_(border);
        // Copy data referred by referrence.
        int the_only_component_index = 0;
        reference_to_component_to_send(component, all_x_, all_y_, all_z_) =
          reference_to_component_border_data(the_only_component_index,
                                             all_x_, all_y_, all_z_);
      }  // end of for component
    }  // end of for border
  }  // end of BasicSimulationCore::PrepareBordersToSend()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Do FDTD stepping for internal part of grid data.
  ///
  /// WARNING! It is repeated MANY times! 
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
    // set all fields (Ex, Ey, Ez, Hx, Hy, Hz) to be zero.
    for (int component = kEx; component <= kHz; ++component) { 
      // data_(component, all_x_, all_y_, all_z_) = 0;
      //debug
      data_(component, all_x_, all_y_, all_z_) = component;      
    }
    blitz::firstIndex component_position_x;
    blitz::secondIndex component_position_y;
    blitz::thirdIndex component_position_z;
    blitz::Array<double, kDimensions> ex = data_(kEx, all_x_, all_y_, all_z_); // view
    /// @todo1 Remove debug assigment for field ex.
    ex = (component_position_y + component_position_x * 10
          + ex.columns() * process_rank_)
      * (component_position_z*2 - 1) * (-1);
    if (process_rank_ == 0) {
      int z = 0;
      std::cout << data_(kEx, all_x_, all_y_, 0) << std::endl;
    }
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
    all_x_ = blitz::Range::all();
    all_y_ = blitz::Range::all();
    all_z_ = blitz::Range::all();
    halo_width_ = simulation_input_config_.halo_width();
    local_time_step_ = 0;
    number_of_components_to_exchange_ = simulation_input_config_.
      number_of_components_to_exchange();
    components_to_exchange_.resize(number_of_components_to_exchange_);
    components_to_exchange_ = simulation_input_config_.components_to_exchange();
    number_of_grid_data_components_ = simulation_input_config_.
      number_of_grid_data_components();
    process_rank_ = process_rank;
    for (int axis = kAxisX; axis < kDimensions; ++axis)
      subdomain_size_[axis] = subdomain_size[axis];
    total_time_steps_ = simulation_input_config_.total_time_steps();
    // Resize grid data to fit subdomain of current MPI process.
    int max_x = subdomain_size_[kAxisX]-1;
    int max_y = subdomain_size_[kAxisY]-1;
    int max_z = subdomain_size_[kAxisZ]-1;
    data_.resize(number_of_grid_data_components_, max_x + 1, max_y + 1, max_z + 1);
    // blitz::RectDomain<kDimensions> borders_range_[kDimensions*2];
    borders_range_.resize(kDimensions*2, number_of_components_to_exchange_);
    for (int component = 0; component < number_of_components_to_exchange_; ++component) {
      // Describe indexes for each border.
      typedef blitz::RectDomain<1+kDimensions> rd;
      typedef blitz::TinyVector<int,4> vec;
      int c = components_to_exchange_(component);
      borders_range_(static_cast<int>(kBorderLeft), component)
        = rd(vec(c,            0,     0,     0),
             vec(c,halo_width_-1, max_y, max_z));
      borders_range_(static_cast<int>(kBorderRight), component)
        = rd(vec(c,max_x - (halo_width_-1),     0,     0),
             vec(c,max_x                  , max_y, max_z));
      borders_range_(static_cast<int>(kBorderBottom), component)
        = rd(vec(c,    0,             0,     0),
             vec(c,max_x, halo_width_-1, max_z));
      borders_range_(static_cast<int>(kBorderTop), component)
        = rd(vec(c,    0, max_y - (halo_width_-1),     0),
             vec(c,max_x, max_y                  , max_z));
      borders_range_(static_cast<int>(kBorderBack), component)
        = rd(vec(c,    0,     0,             0),
             vec(c,max_x, max_y, halo_width_-1));
      borders_range_(static_cast<int>(kBorderFront), component)
        = rd(vec(c,    0,     0, max_z - (halo_width_-1)),
             vec(c,max_x, max_y, max_z                  ));
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
      }
    }
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
    total_time_steps_ = 1;
    // For most simple case we will need Ex, Ey, Ez, epsilon, Hx, Hy,
    // Hz, mu.
    number_of_grid_data_components_ = 8;
    halo_width_ = 2;
    int components_to_exchange[] = {kEx, kEy, kEz, kHx, kHy, kHz};
    number_of_components_to_exchange_ = sizeof(components_to_exchange) / sizeof(int);
    components_to_exchange_.resize(number_of_components_to_exchange_);
    for (int i = 0; i < number_of_components_to_exchange_; ++i)
      components_to_exchange_(i) = components_to_exchange[i];
    /// @brief Nuber of simulation core data components
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
    int64_t length_x = 5, length_y = 9, length_z = 2;
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
