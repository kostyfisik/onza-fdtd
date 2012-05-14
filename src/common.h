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
  /// @brief Error codes
  ///
  /// Error codes used with onza.
  enum Errors {kDone = 0,  ///< no error
               kErrorUsingInputConfigTooEarly  ///< Using InputConfig too early
  };
  /// @brief Status of reading input config.
  ///
  /// For use with onza::SimulationInputConfig::status_
  enum InputConfig {kInputConfigUndefined = 0, kInputConfigAllDone = 10000};
  /// @brief Names for Cartesian axes.
  ///
  /// Use them to access dimesion relative arrays.
  enum Axis {kAxisX = 0, kAxisY, kAxisZ};
  /// @brief Position of halo borders to exchange relative to current process.
  ///
  /// Order in enum is used by GetOppositeBorder()
  enum BorderPosition {kBorderTop = 0, kBorderLeft,  kBorderFront,
                       kBorderBottom,  kBorderRight, kBorderBack};
  /// @brief Returns border opposite to input border
  ///
  /// Calculation of opposite border is based on definition of #BorderPosition.
  /// @param[in] input_border
  /// @return Border opposite to input_border
  inline BorderPosition GetOppositeBorder(BorderPosition input_border) {
    return static_cast<BorderPosition>((input_border+3)%6);
  }
}  // end of namespace onza
#endif  // SRC_COMMON_H_
