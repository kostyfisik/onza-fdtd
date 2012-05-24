#ifndef SRC_SIMULATION_CORE_BASIC_FDTD_H_
#define SRC_SIMULATION_CORE_BASIC_FDTD_H_
///
/// @file   basic-fdtd.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed May  2 12:54:23 2012
///
/// @brief  Simulation data and stepping algorithm
#include <blitz/array.h>
#include "../common.h"
namespace onza {
  // *********************************************************************** //
  /// @brief Input configuration of grid properties.
  class GridInputConfig {
   public:
    /// @brief Default constructor.
    ///
    /// Define zero grid with no elements.
    GridInputConfig() {set_total_grid_length(0);}
    /// @brief Accesor.
    /// @param[in] axis_name Input axis.
    /// @return Size of grid along input axis.
    inline int64_t get_total_grid_length(Axis axis_name) {
      if (axis_name > kDimensions - 1) return 0;
      return total_grid_length_[static_cast<int>(axis_name)];
    };
    /// @brief Mutator.
    /// @param[in] length_x, length_y, length_z Length of grid in
    /// corresponding direction.
    /// @return 0 is OK.
    int set_total_grid_length(int64_t length_x, int64_t length_y = 0,
                              int64_t length_z = 0) {
      if (kDimensions > 0) total_grid_length_[kAxisX] = length_x;
      if (kDimensions > 1) total_grid_length_[kAxisY] = length_y;
      if (kDimensions > 2) total_grid_length_[kAxisZ] = length_z;
      return 0;
    };  // end of set_total_grid_length

   private:
    /// @brief Grid length in each direction.
    ///
    /// Each value is number of grid vertices in corresponding direction.
    int64_t total_grid_length_[kDimensions];
  };  // end of class GridInputConfig
  // *********************************************************************** //
  /// @brief Parsing all input parameters into one object.
  class SimulationInputConfig {
   public:
    /// @brief Constructor. Initially status of reading config is unread.
    SimulationInputConfig():status_(static_cast<int>(kInputConfigUndefined)) {}
    /// @brief Read top level config file
    int ReadConfig();
    /// @brief Grid properties from ReadConfig().
    GridInputConfig grid_input_config_;  // public member to allow direct
                                         // access to GridInputConfig parameters
    /// @brief Accesor to status of simulation input config.
    int status() {return status_;}
    /// @brief Accesor
    int pml_width() {return pml_width_;}
    /// @brief Accesor
    double pml_computational_ratio() {return pml_computational_ratio_;}
    /// @brief Accesor
    int boundary_condition(BorderPosition requested_border) {
      return boundary_condition_[requested_border];}
    /// @brief Accesor
    int halo_width() {return halo_width_;}
    /// @brief Accesor
    int core_data_components() {return core_data_components_;}
    
   private:
    /// @brief Nuber of simulation core data components
    ///
    /// @todo3 Set it automatically from stepping algorithm
    /// description.
    int core_data_components_;
    /// @brief Array of boundary conditions
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions. Values are from
    /// #BoundaryCondition
    int boundary_condition_[6];
    /// @brief Width of halo to exchange.
    int halo_width_;
    /// @brief PML width for model boundary
    ///
    /// Number of grid nodes in PML.
    int pml_width_;
    /// @brief PML nodes computational complexity
    ///
    /// Each node in PML needs more computations compared to regular nodes.
    /// Parameters value is ratio of computational loads from PML node
    /// to regular node. Should be geather than 1.0.
    double pml_computational_ratio_;
    /// @brief Status of reading config with ReadConfig()
    ///
    /// Should be value from #InputConfig
    int status_;
  };  // end of class SimulationInputConfig
  // *********************************************************************** //
  /// @brief Basic class for FDTD simulation.
  ///
  /// Contains all computational domain data, methods to update to next
  /// time step, methods to publish own border and to use foreign borders.
  class BasicSimulationCore {
   public:
    int Init(const int64_t subdomain_size[]);
    /// @brief Parsing all input parameters into one object.
    SimulationInputConfig simulation_input_config_;
   private:
    /// @brief Data components.
    ///
    /// First dimension enumerates data component. Each data component
    /// is a kDimensions array. #DataComponents enum can be used to
    /// access components values, e.g. data_(kEx, x, y, z) or
    /// data_(kEps, x, y, z);
    blitz::Array<double, 1+kDimensions> data_;
  };  // end of class BasicSimulationCore
}  // end of namespace onza
#endif  // SRC_SIMULATION_CORE_BASIC_FDTD_H_
