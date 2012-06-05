#ifndef SRC_COMMON_H_
#define SRC_COMMON_H_
///
/// @file   common.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @date   Wed May  2 13:52:19 2012
///
/// @brief  Global variables, enumerations, macros, types, etc..
#include <stdint.h>  // Including stdint.h to use int64_t.
namespace onza {
  // ********************************************************************** //
  // **********************     Constants           *********************** //
  // ********************************************************************** //
  const int kOutput = 0;
  /// @brief Error codes
  ///
  /// Error codes used with onza.
  enum Errors {
    /// no error
    kDone = 0,
    /// After splitting model to Cartesian grid current procces was
    /// not found in it.
    kErrorProcessNotInGrid,
    /// onza::SimulationInputConfig forced to set total pml width in
    /// some dimension bigger than size of this dimension.
    kErrorTooWidePml,
    /// wrong index while checking with
    /// HaloExchangeProcess::CheckSubdomainIndexes()
    kErrorWrongIndexDifference,
    /// found uninitiated simulation core when
    /// HaloExchangeProcess::RunSimulation()
    kErrorUninitiatedSimulationCore,
    /// case of too many processes during decomposition
    kErrorSubdomainSizeLessThanHaloWidth,
    /// exchange buffer was initiated to be not contiguous with
    /// BasicSimulationCore::Init() or was detected to be not
    /// contiguous with HaloToExchange::Init().
    kErrorExchangeBufferIsNotContiguous,
    /// Buffers to send and recieve  halo has different sizes.
    /// Checked with HaloToExchange::Init().
    kErrorSendAndReceiveBuffersHasDifferentSizes,
    /// Should use some FDTD algorithm.
    kErrorWrongAlgorithm,
    /// FDTD algorithm`s should have some subdomain halo.
    kErrorWrongHaloWidth,
    /// FDTD algorithm`s time depth should be >= 1!
    kErrorWrongTimeDepth,
    /// Using InputConfig too early
    kErrorUsingInputConfigTooEarly
  };
  /// @brief Simulation status
  enum SimulationStatus {kSimulationStatusFinished = 0,
                         kSimulationStatusRunning,
                         kSimulationStatusInitiated,
                         kSimulationStatusUndefined
                         };
  /// @brief Algorithm selection
  enum Algorithm { kAlgorithmSimpleX1D=0, kAlgorithmSimpleY1D,
                   kAlgorithmSimpleZ1D,
                   kAlgorithmSimpleTMz2D,
                   kAlgorithmSimple3D};
  /// @brief Simulation core and halo border predefined names for data
  /// components
  ///
  /// 
  enum DataComponents {kEz=0, kHy, kInvEps, kSrcEz,// 1D axis x
                       kCeze, kCezh, kChyh, kChye,
                       kHx, kChxh, kChxe,          // 2D TM axis z
                       kEx, kCexe, kCexh,
                       kEy, kCeye, kCeyh,
                       kHz, kChzh, kChze};  // 3D FDTD
  /// @brief Status of reading input config.
  ///
  /// For use with onza::SimulationInputConfig::status_
  enum InputConfig {kInputConfigUndefined = 0,
                    kInputConfigErrorWidePml,
                    kInputConfigAllDone = 10000};
  /// @brief Type of boundary conditions for simulation model.
  ///
  /// Currently only PML condition is valid
  enum BoundaryCondition {
    /// Uninitialized boundary condition.
    kBoudaryConditionUndefined = 0,
    /// Boundary is periodical.
    kBoudaryConditionPeriodical,
    /// Boundary is reduced.
    kBoudaryConditionReduced,
    /// There is PML region near this boundary inside current domain.
    kBoudaryConditionPML
  };  // end of enum BoundaryCondition
  /// @brief Total number of Cartesian axes.
  ///
  /// Should allways be 3!!! To define 1D or 2D model set
  /// corresponding grid length(s) to 1.
  const int kDimensions = 3;
  /// @brief Names for Cartesian axes.
  ///
  /// Use them to access dimesion relative arrays.
  enum Axis {kAxisX = 0, kAxisY, kAxisZ};
  /// @brief Position of halo borders to exchange relative to current process.
  ///
  /// Order in enum is used by GetOppositeBorder()
  /// Simple math (due to order in enums): border_axis = border % kDimensions;
  enum BorderPosition {kBorderLeft = 0, kBorderBottom, kBorderBack,
                       kBorderRight, kBorderTop, kBorderFront};
  /// @brief Some values for tags for MPI messages.
  enum MpiTag {kMpiTagCheckIndex = 1};
  // ********************************************************************** //
  // **********************     Global inline functions    **************** //
  // ********************************************************************** //
  /// @brief Returns border opposite to input border
  ///
  /// Calculation of opposite border is based on definition of #BorderPosition.
  /// @param[in] input_border
  /// @return Border opposite to input_border
  inline BorderPosition GetOppositeBorder(BorderPosition input_border) {
    return static_cast<BorderPosition>((input_border+3)%6);
  }
  /// @breif Mutiply components of vector
  ///
  /// Used to get product of any vector type, e.g.
  /// @code
  /// int array_of_int[] = {1, 2, 3};
  /// // 6
  /// printf("Product = %i\n", onza::MultiplyComponents(array_of_int));
  /// double array_of_double[] = { 0.2, 0.3};
  /// // 0.06
  /// printf("Product = %g\n", onza::MultiplyComponents(array_of_double));
  /// @endcode
  /// @param[in] Some vector
  /// @return Product of vector components.
  /// @todo3 check speed of template realization of MultiplyComponents function.
  /// as far as it is belived it has very small performance overhead
  /// (with compiler optimizing flag -O2 or -O3).
  template <class VectorType, int dimensions> inline
      VectorType MultiplyComponents(VectorType(& vector)[dimensions]) {
    if (dimensions == 3) return  vector[0]*vector[1]*vector[2];
    if (dimensions == 2) return  vector[0]*vector[1];
    if (dimensions == 1) return  vector[0];
    if (dimensions > 3) {
      VectorType product = vector[0];
      for (int i = 1; i < dimensions; ++i) {
        product *= vector[i];
      }  // end of for
      return product;
    }  // end of if dimensions
  }  // end of template MultiplyComponents
  /// @breif Sum up components of vector
  ///
  /// Used to get sum of any vector type components.
  /// @see MultiplyComponents()
  /// @param[in] Some vector
  /// @return Sum of vector components.
  /// @todo check speed of template realization of SumUpComponents function.
  /// as far as it is belived it has very small performance overhead
  /// (with compiler optimizing flag -O2 or -O3).
  template <class VectorType, int dimensions> inline
      VectorType SumUpComponents(VectorType(& vector)[dimensions]) {
    if (dimensions == 3) return  vector[0] + vector[1] + vector[2];
    if (dimensions == 2) return  vector[0] + vector[1];
    if (dimensions == 1) return  vector[0];
    if (dimensions > 3) {
      VectorType product = vector[0];
      for (int i = 1; i < dimensions; ++i) {
        product += vector[i];
      }  // end of for
      return product;
    }  // end of if dimensions
  }  // end of template SumUpComponents
  template <class T> inline T pow2(const T value2pow) {return value2pow*value2pow;}
}  // end of namespace onza
#endif  // SRC_COMMON_H_
