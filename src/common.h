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
  /// @brief Error codes
  ///
  /// Error codes used with onza.
  enum Errors {kDone = 0,  /// < no error
               kErrorProcessNotInGrid,  /// < After splitting model to Cartesian
               /// current process was not used.
               kErrorUsingInputConfigTooEarly  /// < Using InputConfig too early
  };
  /// @brief Status of reading input config.
  ///
  /// For use with onza::SimulationInputConfig::status_
  enum InputConfig {kInputConfigUndefined = 0, kInputConfigAllDone = 10000};
  /// @brief Type of boundary conditions for simulation model.
  ///
  /// Currently only PML condition is valid, 
  enum BoundaryCondition {kBoudaryConditionUndefined = 0,
                          kBoudaryConditionHalo,
                          kBoudaryConditionPML};
  /// @brief Total number of Cartesian axes.
  ///
  /// It is tested for 3D models, sometimes 1D and 2D are also tested.
  /// Valid values are only 1, 2 and 3 (dimensions).
  const int kDimensions = 3;
  /// @brief Names for Cartesian axes.
  ///
  /// Use them to access dimesion relative arrays.
  enum Axis {kAxisX = 0, kAxisY, kAxisZ};
  /// @brief Position of halo borders to exchange relative to current process.
  ///
  /// Order in enum is used by GetOppositeBorder()
  enum BorderPosition {kBorderRight = 0, kBorderTop, kBorderFront,
                       kBorderLeft, kBorderBottom, kBorderBack};
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
}  // end of namespace onza
#endif  // SRC_COMMON_H_
