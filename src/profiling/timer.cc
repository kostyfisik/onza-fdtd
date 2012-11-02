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
#include "timer.h"
namespace onza {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
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
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
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
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Single MPI_Wtime() call for sequental timers.
  int Timer::StopStart(const char *stop_key_input, const char *start_key_input) {
    std::string stop_key(stop_key_input), start_key(start_key_input);
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::Print(const char *key_input) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknow timer!\n");      
      return 1; //kErrorProfilingTimerSecondStart;
    }
    printf("\t%04.2f s : %s\n", total_time_[key], key.c_str());
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::PrintRelative(const char *key_input, double base) {
    std::string key(key_input);
    if (!total_time_.count(key)) {
      printf("Error! Printing unknow timer!\n");      
      return 1; //kError
    }
    if (base == 0) {
      printf("Error! Printing timer with zero base!\n");      
      return 1; //kError
    }    
    double time = total_time_[key];
    printf("\t%04.1f%%\t%3.2f s\t%s\n", time/base*100.0, time,  key.c_str());
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::PrintAll() {
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::PrintAllRelative(double base) {
    if (base == 0) {
      printf("Error! Printing timer with zero base!\n");      
      return 1; //kError
    }    
    std::map<std::string, double>::iterator it;
    // Sorted by value, sorting keys.
    std::map<std::string, std::string> sorted_keys;
    for (it=total_time_.begin() ; it != total_time_.end(); ++it) {
      std::string sorting_key = to_string((*it).second/base*100.0)
        + (*it).first;
      sorted_keys[sorting_key] = (*it).first;
    }  // end of for sorting keys
    // Sorted by value, output.
    std::map<std::string, std::string>::reverse_iterator it2;
    for (it2=sorted_keys.rbegin() ; it2 != sorted_keys.rend(); ++it2)
      PrintRelative((*it2).second.c_str(),base);
    // // Simple output, sorted by key name.
    // printf("Sorted by key name.\n");
    // for (it=total_time_.begin() ; it != total_time_.end(); ++it)
    //   PrintRelative((*it).first.c_str(),base);
    // return kDone;
  }  // end of int Timer::PrintAllRelative()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  double Timer::GetKeyTotalTime(const char *key_input) {
    std::string key(key_input);
    return total_time_[key];
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  double Timer::GetAllTotalTime() {
    std::map<std::string, double>::iterator it;
    double all_total_time = 0;
    for (it=total_time_.begin() ; it != total_time_.end(); ++it)
      all_total_time += (*it).second;
    return all_total_time;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::Reset(const char *key_input) {
    std::string key(key_input);
    return kDone;
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief
  int Timer::ResetAll() {
    return kDone;
  }
}
