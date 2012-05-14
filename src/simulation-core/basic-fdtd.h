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
    GridInputConfig() {set_total_grid_length(0, 0, 0);}
    /// @brief Accesor.
    /// @param[in] axis_name Input axis.
    /// @return Size of grid along input axis.
    inline int64_t get_total_grid_length(Axis axis_name) {
      return total_grid_length_[static_cast<int>(axis_name)];
    };
    /// @brief Mutator.
    /// @param[in] length_x, length_y, length_z Length of grid in
    /// corresponding direction.
    /// @return 0 is OK.
    int set_total_grid_length(int64_t length_x, int64_t length_y,
                              int64_t length_z) {
      total_grid_length_[0] = length_x;
      total_grid_length_[1] = length_y;
      total_grid_length_[2] = length_z;
      return 0;
    };  // end of set_total_grid_length

   private:
    /// @brief Grid length in each direction.
    int64_t total_grid_length_[3];
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
   private:
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
    int Init();
    /// @brief Parsing all input parameters into one object.
    SimulationInputConfig simulation_input_config_;
   private:
    /// @brief Some array for construction tests only.
    blitz::Array<double, 3> test_array_;
  };  // end of class BasicSimulationCore
}  // end of namespace onza
#endif  // SRC_SIMULATION_CORE_BASIC_FDTD_H_
