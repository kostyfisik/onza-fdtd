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
#include <map>
#include <string>
#include "../common.h"
namespace onza {
  // *********************************************************************** //
  /// @brief Input configuration of grid properties.
  class GridInputConfig {
   public:
    /// @brief Default constructor.
    ///
    /// Define minimal grid with single element.
    GridInputConfig() {set_total_grid_length(1, 1, 1);}
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
    int set_total_grid_length(int64_t length_x, int64_t length_y,
                              int64_t length_z) {
      if (length_x < 1) return kErrorSettingWrongGridSize;
        total_grid_length_[kAxisX] = length_x;
      if (length_y < 1) return kErrorSettingWrongGridSize;
      total_grid_length_[kAxisY] = length_y;
      if (length_z < 1) return kErrorSettingWrongGridSize;
      total_grid_length_[kAxisZ] = length_z;
      return kDone;
    };  // end of set_total_grid_length

   private:
    /// @brief Grid length in each direction.
    ///
    /// Each value is number of grid vertices in corresponding direction.
    /// Should allways be positive (>=1), as far as it is used both
    /// to evaluate size of simualtion domain and describe virtual
    /// topology of MPI processe for halo exchange.
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
    int time_depth() {return time_depth_;}
    /// @brief Accesor
    int number_of_grid_data_components()
      {return number_of_grid_data_components_;}
    /// @brief Accesor
    int number_of_components_to_exchange()
      {return number_of_components_to_exchange_;}
    /// @brief Accesor
    blitz::Array<int, 1> components_to_exchange()
      {return components_to_exchange_;}
    /// @brief Accesor
    double snapshot_interval() {return snapshot_interval_;}
    /// @brief Accesor
    int64_t total_time_steps() {return total_time_steps_;}
    /// @brief Accesor
    int algorithm() {return algorithm_;}
    int set_config_file_name(std::string config_file_name) {
      config_file_name_ = config_file_name;
      return kDone;
    };

   private:
    /// @brief Auto set boundary conditions for reduced dimensions
    int AutoSetReducedBoundaryConditions();
    /// @brief Check total PML width to be less than domain width
    int CheckTotalPMLWidth();
    /// @brief Parse config file keys.
    int SetConfigFileMap();
    /// @brief Algorithm selection (from #Algorithm)
    int algorithm_;
    /// @brief Array of boundary conditions
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions. Values are from
    /// #BoundaryCondition
    int boundary_condition_[6];
    /// @brief List of components to exchange in halo.
    blitz::Array<int, 1> components_to_exchange_;
    /// @brief Configuration file name.
    std::string config_file_name_;
    /// @brief Maped values from configuration file.
    std::map<std::string, std::string> config_file_map_;
    /// @brief Number of simulation core data components to exchange
    /// in halo
    int number_of_components_to_exchange_;
    /// @brief Nuber of simulation core data components
    ///
    /// @todo3 Set it automatically from stepping algorithm
    /// description.
    int number_of_grid_data_components_;
    /// @brief Depth of FDTD in time.
    ///
    /// Simple FDTD has depth equal to 2 - to get data for next step
    /// it needs data from previous step.
    int time_depth_;
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
    /// @brief Total time steps in simulation.
    int64_t total_time_steps_;
    /// @brief Status of reading config with ReadConfig()
    ///
    /// Should be value from #InputConfig
    int status_;
    /// @breif snapshot interval
    double snapshot_interval_;
  };  // end of class SimulationInputConfig
  // *********************************************************************** //
  /// @brief Basic class for FDTD simulation.
  ///
  /// Contains all computational domain data, methods to update to next
  /// time step, methods to publish own border and to use foreign borders.
  class BasicSimulationCore {
   public:
    /// @name Borders.
    /// @warning This members are direclty accessed with
    /// onza::HaloToExchange by pointer arithmetics. Be sure to
    /// synchronize reads and writes with it. It is intended to be
    /// done with HaloExchangeProcess::RunSimulation().
    // @{
    /// @bief Current process borders to be send to neighbours.
    ///
    /// First dim - border name (from kBorderLeft to kBorderFront)
    /// Second dim - grid data component
    /// Last three dims - kAxisX, kAxisY, kAxisZ.
    /// @warning This member is direclty accessed with
    /// onza::HaloToExchange by pointer arithmetics. Be sure to
    /// synchronize reads and writes with it. It is intended to be
    /// done with HaloExchangeProcess::RunSimulation().
    blitz::Array<blitz::Array<double, 1+kDimensions>, 1> borders_to_send_;
    /// @brief Received borders
    ///
    /// Borders, received from neighbours.
    /// First dim - border name (from kBorderLeft to kBorderFront)
    /// Second dim - grid data component
    /// Last three dims - kAxisX, kAxisY, kAxisZ.
    /// @warning This member is direclty accessed with
    /// onza::HaloToExchange by pointer arithmetics. Be sure to
    /// synchronize reads and writes with it. It is intended to be
    /// done with HaloExchangeProcess::RunSimulation().
    blitz::Array<blitz::Array<double, 1+kDimensions>, 1> received_borders_;
    // @}
    /// @brief Parsing all input parameters into one object.
    SimulationInputConfig simulation_input_config_;
    /// @breif Initialize member data.
    int Init(const int64_t subdomain_size[],
             int process_rank, int neighbours_ranks[],
             int my_coords[],
             int64_t subdomain_start_index[],
             int64_t subdomain_finish_index[]);
    /// @brief Initialize grid data for all components
    int SetGridData();
    /// @brief Prepare borders to send
    void PrepareBordersToSend();
    /// @brief Snapshot component from data_ to fiel.
    void Snapshot();
    /// @brief Prepare source component in data_.
    void PrepareSource();
    /// @brief Cycle snapshots before new timestep.
    void CycleSnapshots();
    /// @brief Do FDTD stepping for border part of grid data.
    void DoBorderStep();
    /// @brief Do FDTD stepping for internal part of grid data.
    void DoStep();
    /// @brief
    int StepTime();
    /// @brief Accesor.
    int status() {return status_;}
    int64_t total_time_steps() {return total_time_steps_;}

   private:
    /// @brief Algorithm selection (from #Algorithm)
    int algorithm_;
    /// @brief FDTD algorithm for 3D case.
    void AlgorithmSimple3D(blitz::Range x,
                           blitz::Range y,
                           blitz::Range z);
    /// @brief FDTD algorithm for TMz 2D case.
    void AlgorithmSimpleTMz2D(blitz::Range x,
                              blitz::Range y,
                              blitz::Range z);
    /// @brief FDTD algorithm for 1D case (length in kAxisX dimenstion).
    void AlgorithmSimpleX1D(blitz::Range x,
                            blitz::Range y,
                            blitz::Range z);
    /// @brief FDTD algorithm for 1D case (length in kAxisY dimenstion).
    void AlgorithmSimpleY1D(blitz::Range x,
                            blitz::Range y,
                            blitz::Range z);
    /// @brief FDTD algorithm for 1D case (length in kAxisZ dimenstion).
    void AlgorithmSimpleZ1D(blitz::Range x,
                            blitz::Range y,
                            blitz::Range z);
    /// @brief Predefined ranges for all grid points in x, y and z dimensions.
    blitz::Range all_x_, all_y_, all_z_;
    /// @brief Predefined ranges for data grid points inside border.
    blitz::Range data_border_range_[kDimensions*2];
    /// @brief Depth of FDTD in time.
    ///
    /// Simple FDTD has depth equal to 2 - to get data for next step
    /// it needs data from previous step.
    int time_depth_;
    /// @brief Width of halo to exchange.
    int halo_width_;
    /// @brief Predefined ranges for inner grid points in x, y and z dimensions.
    ///
    /// Inner - means that next step for inner point could be
    /// calculated without using borders information.
    blitz::Range inner_x_, inner_y_, inner_z_;
    /// @brief Current time step for current MPI process.
    int64_t local_time_step_;
    /// @brief This process coords in Cartesian topology of MPI processes.
    int my_coords_[kDimensions];
    /// @brief My neighbours ranks
    ///
    /// Due to order in enum #BorderPosition using max(kDimensions)*2
    /// size of array for all kDimensions.
    /// @see HaloExchangeProcess::neighbours_ranks_
    int neighbours_ranks_[6];
    /// @brief Number of simulation core data components to exchange
    /// in halo
    int number_of_components_to_exchange_;
    /// @brief Number of simulation core data components
    int number_of_grid_data_components_;
    /// @brief Fast access (without MPI call) to process rank of
    /// containing object. Should be set in init()
    int process_rank_;
    /// @brief Pointer to FDTD algorithm used for stepping.
    ///
    /// FDTD algorithm should be selected during Init()
    void (BasicSimulationCore::*RunAlgorithm)(blitz::Range x,
                                              blitz::Range y,
                                              blitz::Range z);
    /// @breif snapshot interval
    double snapshot_interval_;
    int snapshot_frame_;
    /// @brief Simulation current.
    int status_;
    /// @brief Current process subdomain size in each dimension.
    ///
    /// size == number of vertices
    int64_t subdomain_size_[kDimensions];
    /// @breif Index of global model pointing to the beginnig of the
    /// subdomain.
    ///
    /// 0 index of subdomain should refer to same part of model
    /// as subdomain_start_index_
    int64_t subdomain_start_index_[kDimensions];
    /// @breif Index of global model pointing to the end of the
    /// subdomain.
    ///
    /// @breif subdomain_size_-1 index of subdomain should refer to
    /// same part of model as subdomain_finish_index_
    int64_t subdomain_finish_index_[kDimensions];
    /// @brief Total time steps in simulation.
    int64_t total_time_steps_;
    /// @brief List of components to exchange in halo.
    blitz::Array<int, 1> components_to_exchange_;
    /// @brief Grid data components.
    ///
    /// First dimension enumerates data component. Each data component
    /// is a kDimensions array. #DataComponents enum can be used to
    /// access components values, e.g. data_(kEx, x, y, z) or
    /// data_(kEps, x, y, z);
    ///
    /// @todo1 Change from double to custom (for ease of switching
    /// from double to single precision). This will affect cache usage
    /// (during stepping algorithm) and network utilization (during
    /// halo exchange process).
    ///
    /// @todo1 Change access pattern from data_(kEx, x, y, z) to
    /// data_(x, y, z, kEx). This can icrease cache usage, especially
    /// for 3D case
    blitz::Array<double, 1+kDimensions> data_;
    /// @brief Grid data components snapshots at different timesteps.
    ///
    /// data_snapshot_(time_depth_ - 2) is a reference to data_ and an only
    /// data_snapshot_ exchanging its borders.
    /// data_snapshot_(time_depth_ - 1) is reserved for the results,
    /// calculated at current step.
    blitz::Array<blitz::Array<double, 1+kDimensions>, 1> data_snapshot_;
    /// @brief Borders ranges inside grid.
    /// First dim - border name (from kBorderLeft to kBorderFront)
    /// Second dim - grid data component.
    /// e.g borders_to_send_range_(kBorderLeft, 2) is range of indexes in
    /// data_ for border kBorderLeft and component
    /// components_to_exchange_(2).
    blitz::Array<blitz::RectDomain<1+kDimensions>, 2> borders_to_send_range_,
        received_borders_range_;
  };  // end of class BasicSimulationCore
}  // end of namespace onza
#endif  // SRC_SIMULATION_CORE_BASIC_FDTD_H_
