///
/// @file   timer.cc
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
  /// @brief Set a start mark.
  /// Sets a start mark for the input element(e.g. any function name
  /// to do profiling).
  int Timer::Start(const char *key_input) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      total_time_[key] = 0;
      is_running_[key] = kTimerOff;
    }
    if (is_running_[key] != kTimerOff) {
      printf("Error! Starting timer second time!\n");
      return kErrorProfilingTimerSecondStart;
    }
    is_running_[key] = kTimerOn;
    start_mark_[key] = MPI_Wtime();
    return kDone;
  } // end of Timer::Start
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Set a stop mark and measure time.
  /// Sets a stop mark for the input element(e.g. any function) and
  /// adds total time (stop mark - start mark) spent on element's processing.
  int Timer::Stop(const char *key_input) {
    std::string key(key_input);
    stop_mark_[key] = MPI_Wtime();
    if (is_running_[key] != kTimerOn || !total_time_.count(key)) {
      printf("Error! Wrong stopping timer!\n");
      return kErrorProfilingTimerWrongStop;
    }
    is_running_[key] = kTimerOff;
    total_time_[key] += stop_mark_[key] - start_mark_[key];
    return kDone;
  } // end of Timer::Stop
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Stop timer for the first input parameter and start for the second.
  /// Stops timer for the first input parameter (...stop_key_input) and
  /// starts it for the second one (...start_key_input).
  int Timer::StopStart(const char *stop_key_input,
                       const char *start_key_input) {
    std::string stop_key(stop_key_input), start_key(start_key_input);
    if (is_running_[stop_key] != kTimerOn) {
      printf("Error! Attempting to stop the timer that is off!\n");
      printf("Check input arguments to be (stop_key, start_key)\n");
      printf("or consider using Timer::Stop to stop timer properly.\n");
    }
    else {
      stop_mark_[stop_key] = MPI_Wtime();
      total_time_[stop_key] += stop_mark_[stop_key] - start_mark_[stop_key];
      is_running_[stop_key] = kTimerOff;
      if (is_running_[start_key] != kTimerOff) {
        printf("Error! Attempting to start the timer that is on!\n");
        printf("Check input arguments to be (stop_key, start_key)\n");
        printf("or consider using Timer::Start to start timer properly.\n");   
      }
      else {
        start_mark_[start_key] = stop_mark_[stop_key];
        is_running_[start_key] = kTimerOn;
      }
    }
    return kDone;
  } // end of Timer::StopStart
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print total time spent on processing element.
  /// Prints calculated total time for input parameter.
  int Timer::Print(const char *key_input) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknown timer!\n");
      return 1;  // kErrorProfilingTimerSecondStart;
    }
    printf("\t%04.2f s : %s\n", total_time_[key], key.c_str());
    return kDone;
  }  // end of Timer::Print
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Print total time spent on processing element as ratio to a base.
  /// Prints total time spent on processing element as a ratio to a base.
  int Timer::PrintRelative(const char *key_input, double base) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknown timer!\n");
      return 1;  // kError
    }
    if (base == 0) {
      printf("Error! Printing timer with zero base!\n");
      return 1;  // kError
    }
    double time = total_time_[key];
    printf("\t%04.1f%%\t%3.2f s\t%s\n", time/base*100.0, time,  key.c_str());
    return kDone;
  } // end of Timer::PrintRelative
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief For each element print its key(name) and value(time, spent on
  /// processing).
  /// Prints total time spent on processing each element.
  int Timer::PrintAll() {
     std::map<std::string, double>::iterator it;
    // Sorted by value, sorting keys.
    std::map<std::string, std::string> sorted_keys;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      std::string sorting_key = to_string((*it).second * 1.0)
        + (*it).first;
      sorted_keys[sorting_key] = (*it).first;
    }  // end of for sorting keys
    // Sorted by value, output.
    std::map<std::string, std::string>::reverse_iterator it2;
    for (it2 = sorted_keys.rbegin(); it2 != sorted_keys.rend(); ++it2)
      Print((*it2).second.c_str());
    return kDone;
  } // end of Timer::PrintAll
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief For each element print its key(name) and value(time, spent on
  /// processing the element, divided by base) sorted by value.
  /// Prints total time spent on processing each element in percentage to base.
  int Timer::PrintAllRelative(double base) {
    if (base == 0) {
      printf("Error! Printing timer with zero base!\n");
      return kErrorZeroBase;  // kError
    }
    std::map<std::string, double>::iterator it;
    // Sorted by value, sorting keys.
    std::map<std::string, std::string> sorted_keys;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      std::string sorting_key = to_string((*it).second/base*100.0)
        + (*it).first;
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
  /// @brief Print all elements' values sorted by key.
  /// Prints all elements' values sorted by key.
  int Timer::PrintByKey () {
    std::map<std::string, double>::iterator it;
    for (it = total_time_.begin() ; it != total_time_.end(); ++it)
      Print((*it).first.c_str());
    return kDone;
  } // end of Timer::PrintByKey ()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Get total time spent on processing the element.
  /// Gets total time spent on processing the element.
  double Timer::GetKeyTotalTime(const char *key_input) {
    std::string key(key_input);
    return total_time_[key];
  } // end of Timer::GetKeyTotalTime
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Get total time spent on processing all elements.
  /// Gets total time spent on processing all elements.
  double Timer::GetAllTotalTime() {
    std::map<std::string, double>::iterator it;
    double all_total_time = 0;
    for (it = total_time_.begin(); it != total_time_.end(); ++it)
      all_total_time += (*it).second;
    return all_total_time;
  } // end of Timer::GetAllTotalTime
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Reset timer for selected element
  int Timer::Reset(const char *key_input) {
    std::string key(key_input);
    if (is_running_[key] != kTimerOn) {
      printf("Attempting to reset the timer that is off!\n");
    }
    else {
      is_running_[key] = kTimerOff;
      Start(key.c_str());
    }
    return kDone;
  } // end of Timer::Reset
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Reset timer for all elements
  int Timer::ResetAll() {
    std::map<std::string, double>::iterator it;
    for (it = total_time_.begin(); it != total_time_.end(); ++it) {
      if (is_running_[(*it).first] != kTimerOn) {
        printf("Attempting to reset the timer that is off!\n");
        printf("Key = %s\n",(*it).first.c_str());
      }
      else { Start((*it).first.c_str()); }
    }
    return kDone;
  } // end of Timer::ResetAll
}
