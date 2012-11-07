#ifndef SRC_PROFILING_TIMER_H_
#define SRC_PROFILING_TIMER_H_
///
/// @file   timer.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of Onza FDTD.
///
/// Onza FDTD is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// Onza FDTD is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Onza FDTD.  If not, see <http://www.gnu.org/licenses/>.
/// @date   Thu Nov  1 13:36:49 2012
///
/// @brief  Number of tools for profiling.
#include <map>
#include <string>
namespace onza {
  /// @brief MPI_Wtime based timer for profiling.
  class Timer {
   public:
    /// @brief Init (if needed) timer and set a start time mark.
    int Start(const char *key_input);
    /// @brief Set a stop time mark and measure time.
    int Stop(const char *key_input);
    /// @brief Sequentially stop and start timers.
    int StopStart(const char *stop_key_input, const char *start_key_input);
    /// @brief Print total time spent on processing element.
    int Print(const char *key_input);
    /// @brief Print total time spent on processing element as ratio to a base.
    int PrintRelative(const char *key_input, double base);
    /// @brief Print all collected data (ordered).
    int PrintAll();
    /// @brief Print all collected data (ordered, percentage of base time).
    int PrintAllRelative(double base);
    /// @brief Print all timer values sorted by key.
    int PrintAllByKey();
    /// @brief Get total time spent on processing the element.
    double GetKeyTotalTime(const char *key_input);
    /// @brief Get total time spent on processing all elements.
    double GetAllTotalTime();
    /// @brief Preset timer.
    int Preset(const char *key_input, const double total_time);
    /// @brief Reset timer.
    int Reset(const char *key_input);
    /// @brief Reset timer for all elements.
    int ResetAll();

   private:
    /// @brief Timer start time mark.
    std::map<std::string, double> start_mark_;
    /// @brief Timer stop time mark.
    std::map<std::string, double> stop_mark_;
    /// @brief Total time (stop - start).
    std::map<std::string, double> total_time_;
    /// @brief Timer status.
    std::map<std::string, int> is_running_;
    enum TimerStatus {kTimerOff = 0, kTimerOn};
  };
}
#endif  // SRC_PROFILING_TIMER_H_
