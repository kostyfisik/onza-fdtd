///
/// @file   timer.cc
/// @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
/// @copyright 2012 Konstantin Ladutenko
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
/// @date   Thu Nov  1 14:10:45 2012
///
/// @brief MPI_Wtime based timer for profiling
#include <mpi.h>
#include <map>
#include <string>
#include "../common.h"
#include "./timer.h"
namespace onza {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Init (if needed) timer and set a start time mark.
  ///
  /// Sets a start mark for the input element (e.g. any function name
  /// to do profiling).
  int Timer::Start(const char *key_input) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      total_time_[key] = 0;
      is_running_[key] = kTimerOff;
    }  // end of if starting timer for new key
    if (is_running_[key] != kTimerOff) {
      printf("Error! Starting timer %s second time!\n", key.c_str());
      return kErrorProfilingTimerSecondStart;
    }  // end of if starting timer already running
    is_running_[key] = kTimerOn;
    start_mark_[key] = MPI_Wtime();
    return kDone;
  }  // end of Timer::Start()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Set a stop time mark and measure time.
  ///
  /// Sets a stop mark for the input element (e.g. any function) and
  /// adds total time (stop mark - start mark) spent on element's
  /// processing.
  int Timer::Stop(const char *key_input) {
    std::string key(key_input);
    stop_mark_[key] = MPI_Wtime();
    if (is_running_[key] != kTimerOn || !total_time_.count(key)) {
      printf("Error! Stopping undefinded timer %s!\n", key.c_str());
      return kErrorProfilingTimerWrongStop;
    }  // end of if trying to stop already stopped or not existing timer
    is_running_[key] = kTimerOff;
    total_time_[key] += stop_mark_[key] - start_mark_[key];
    return kDone;
  }  // end of Timer::Stop()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Sequentially stop and start timers.
  ///
  /// Stops timer for the first input parameter (...stop_key_input)
  /// and starts it for the second one (...start_key_input).
  /// Logically it is equal to sequential call of Stop() and Start(),
  /// but uses single call of MPI_Wtime().
  int Timer::StopStart(const char *stop_key_input,
                       const char *start_key_input) {
    std::string stop_key(stop_key_input), start_key(start_key_input);
    double time_mark = MPI_Wtime();
    // Stop.
    stop_mark_[stop_key] = time_mark;
    if (is_running_[stop_key] != kTimerOn || !total_time_.count(stop_key)) {
      printf("Error! Stopping undefinded timer %s!\n", stop_key.c_str());
      return kErrorProfilingTimerWrongStop;
    }  // end of if trying to stop already stopped or not existing timer
    is_running_[stop_key] = kTimerOff;
    total_time_[stop_key] += stop_mark_[stop_key] - start_mark_[stop_key];
    // Start.
    if (!total_time_.count(start_key)) {
      total_time_[start_key] = 0;
      is_running_[start_key] = kTimerOff;
    }  // end of if starting timer for new start_key
    if (is_running_[start_key] != kTimerOff) {
      printf("Error! Starting timer %s second time!\n", start_key.c_str());
      return kErrorProfilingTimerSecondStart;
    }  // end of if starting timer already running
    is_running_[start_key] = kTimerOn;
    start_mark_[start_key]  = time_mark;
    return kDone;
  }  // end of Timer::StopStart()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print total time spent on processing element.
  ///
  /// Prints calculated total time for input parameter.
  int Timer::Print(const char *key_input) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknown timer %s!\n", key.c_str());
      return kErrorProfilingTimerPrintUnknownTimer;
    }  // end of trying to print timer data for nonexistent timer
    if (is_running_[key] != kTimerOff) {
      printf("Error! Timer %s is not ready for printing!\n", key.c_str());
      return kErrorProfilingTimerSecondStart;
    }  // end of if timer already running
    printf("\t%13.2f s : %s\n", total_time_[key], key.c_str());
    return kDone;
  }  // end of Timer::Print()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print total time spent on processing element as ratio to a base.
  ///
  /// Prints total time spent on processing element as a ratio to a base.
  int Timer::PrintRelative(const char *key_input, double base) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknown timer %s!\n", key.c_str());
      return kErrorProfilingTimerPrintUnknownTimer;
    }  // end of trying to print timer data for nonexistent timer
    if (is_running_[key] != kTimerOff) {
      printf("Error! Timer %s is not ready for printing!\n", key.c_str());
      return kErrorProfilingTimerSecondStart;
    }  // end of if timer already running
    if (base == 0) {
      printf("Error! Printing timer %s with zero base!\n", key.c_str());
      return kErrorProfilingTimerZeroBase;
    }  // end of if base is zero
    double time = total_time_[key];
    printf("\t%04.1f%%\t%5.2f s : %s\n", time/base*100.0, time,  key.c_str());
    return kDone;
  }  // end of Timer::PrintRelative()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print all collected data (ordered).
  ///
  /// Prints total time spent on processing each element. Sorted by
  /// value, largest values first.
  int Timer::PrintAll() {
    double base = GetAllTotalTime();
    std::map<std::string, double>::iterator it;
    // Sorted by value, sorting keys.
    std::map<std::string, std::string> sorted_keys;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      std::string sorting_key =
        to_string((*it).second / base * 100.0) + (*it).first;
      sorted_keys[sorting_key] = (*it).first;
    }  // end of for sorting keys
    // Sorted by value, output.
    std::map<std::string, std::string>::reverse_iterator it2;
    for (it2 = sorted_keys.rbegin(); it2 != sorted_keys.rend(); ++it2)
      Print((*it2).second.c_str());
    return kDone;
  }  // end of Timer::PrintAll()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print all collected data (ordered, percentage of base time).
  ///
  /// Prints total time spent on processing each element in percentage
  /// of base time. Sorted by value, largest values first.
  int Timer::PrintAllRelative(double base) {
    if (base == 0) {
      printf("Error! Printing timers with zero base!\n");
      return kErrorProfilingTimerZeroBase;
    }  // end of if base is zero
    std::map<std::string, double>::iterator it;
    // Sorting keys by value using map internal layout.
    std::map<std::string, std::string> sorted_keys;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      std::string sorting_key =
        to_string((*it).second / base * 100.0) + (*it).first;
      sorted_keys[sorting_key] = (*it).first;
    }  // end of for sorting keys
    // Sorted by value, output.
    std::map<std::string, std::string>::reverse_iterator it2;
    for (it2 = sorted_keys.rbegin(); it2 != sorted_keys.rend(); ++it2)
      PrintRelative((*it2).second.c_str(), base);
    return kDone;
  }  // end of Timer::PrintAllRelative()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print all timer values sorted by key.
  ///
  /// Prints total time spent on processing each element in
  /// seconds. Sorted by timer key name.
  int Timer::PrintAllByKey() {
    std::map<std::string, double>::iterator it;
    for (it = total_time_.begin() ; it != total_time_.end(); ++it)
      Print((*it).first.c_str());
    return kDone;
  }  // end of Timer::PrintAllByKey()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Get total time spent on processing the element.
  ///
  /// Gets total time spent on processing the element.
  double Timer::GetKeyTotalTime(const char *key_input) {
    std::string key(key_input);
    return total_time_[key];
  }  // end of Timer::GetKeyTotalTime()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Get total time spent on processing all elements.
  ///
  /// Gets total time spent on processing all elements.
  double Timer::GetAllTotalTime() {
    std::map<std::string, double>::iterator it;
    double all_total_time = 0;
    for (it = total_time_.begin(); it != total_time_.end(); ++it)
      all_total_time += (*it).second;
    return all_total_time;
  }  // end of Timer::GetAllTotalTime()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Preset timer.
  ///
  /// Set timer total_time_ to some preset value.
  int Timer::Preset(const char *key_input, const double total_time) {
    std::string key(key_input);
    if (is_running_[key] != kTimerOff) {
      printf("Attempting to preset the timer %s, while it is on!\n",
             key.c_str());
      return kErrorProfilingTimerWrongPreset;
    }  // end of if timer was running
    total_time_[key] = total_time;
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Reset timer.
  ///
  /// Restores timer to it uninitialized status.
  int Timer::Reset(const char *key_input) {
    std::string key(key_input);
    if (is_running_[key] != kTimerOff) {
      printf("Attempting to reset the timer %s, while it is on!\n",
             key.c_str());
      return kErrorProfilingTimerWrongReset;
    }  // end of if timer was running
    start_mark_.erase(key);
    stop_mark_.erase(key);
    total_time_.erase(key);
    is_running_.erase(key);
    return kDone;
  }  // end of Timer::Reset()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Reset timer for all elements.
  ///
  /// Resets timer for all elements
  int Timer::ResetAll() {
    std::map<std::string, double>::iterator it;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      if (is_running_[(*it).first] != kTimerOff) {
        printf("Attempting to reset the timer %s, while it is on!\n",   \
               (*it).first.c_str());
        return kErrorProfilingTimerWrongReset;
      }  // end of if timer was running
    }  // end of for all timers check running status
    start_mark_.clear();
    stop_mark_.clear();
    total_time_.clear();
    is_running_.clear();
    return kDone;
  }  // end of Timer::ResetAll()
}
