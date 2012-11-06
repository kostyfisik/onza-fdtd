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
  /// @brief MPI_Wtime based timer for profiling
  class Timer {
  public:
    int Start(const char *key_input);
    int Stop(const char *key_input);
    int StopStart(const char *stop_key_input, const char *start_key_input);
    int Print(const char *key_input);
    int PrintRelative(const char *key_input, double base);
    int PrintAll();
    int PrintAllRelative(double base);
    double GetKeyTotalTime(const char *key_input);
    double GetAllTotalTime();
    int Reset(const char *key_input);
    int ResetAll();
  private:
    std::map<std::string, double> start_mark_;
    std::map<std::string, double> stop_mark_;
    std::map<std::string, double> total_time_;
    std::map<std::string, int> is_running_;
    enum TimerStatus {kTimerOff = 0, kTimerOn};
  };
}
#endif  // SRC_PROFILING_TIMER_H_
