/// @file   2D_test_range.cc
/// @author Markovich Dmitry <dmmrkovich at gmail (.) com>
/// @copyright 2013 Markovich Dmitry
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
#include <blitz/array.h>
#include <blitz/array/stencilops.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include <ctime>

#define SIZE sz
#define STEPS st

// Setting model parameters
const int size = SIZE;
const int steps = STEPS;
const int length_x = SIZE;
const int length_y = SIZE;
const int repeats = 2;
const int sum_flag = 0, time_flag = 1, print_flag = 0;
const double initial = 0.0, coefficient = 0.1, one = 1.0;
// End Setting model parameters

blitz::Array <double, 3> data0, data1;
const blitz::Range i = blitz::Range::all(), j = blitz::Range::all();

enum DataComponents {kEz = 0, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kSrcEz};
const int components = kSrcEz + 1;

void UpdateFields2DRange(blitz::Range temp_x, blitz::Range temp_y) {
  data1(kHx, temp_x, temp_y) = data0(kChxh, temp_x, temp_y) *
    data0(kHx, temp_x, temp_y) - data0(kChxe, temp_x, temp_y) *
    (data0(kEz, temp_x, temp_y+1) - data0(kEz, temp_x, temp_y));

  data1(kHy, temp_x, temp_y) = data0(kChyh, temp_x, temp_y) *
    data0(kHy, temp_x, temp_y) + data0(kChye, temp_x, temp_y) *
    (data0(kEz, temp_x+1, temp_y) - data0(kEz, temp_x, temp_y));

  data1(kEz, temp_x, temp_y) = data0(kCeze, temp_x, temp_y) *
    data0(kEz, temp_x, temp_y) + data0(kCezh, temp_x, temp_y) *
    ((data0(kHy, temp_x, temp_y) - data0(kHy, temp_x-1, temp_y))
     - (data0(kHx, temp_x, temp_y) - data0(kHx, temp_x, temp_y-1)))
    + data0(kSrcEz, temp_x, temp_y);
}

double UsingRange(int length_x, int length_y, int steps) {
  const blitz::Range inner_i(1, length_x);
  const blitz::Range inner_j(1, length_y);
  double time = 0, start, end;
  for (int l = 0; l < repeats; ++l) {
    data0(kEz, i, j) = initial;
    data0(kHx, i, j) = initial;
    data0(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src = exp(-static_cast<double>(t * t) / (steps * steps));
      data0(static_cast<int>(kSrcEz), i, j) = src;
      UpdateFields2DRange(inner_i, inner_j);
      cycleArrays(data0, data1);
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Range 2D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Range 2D control sum for Ez is "
              << sum(abs(data0(static_cast<int>(kEz), i, j))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data0(static_cast<int>(kHx), i, j);
    std::cout << data0(static_cast<int>(kHy), i, j);
    std::cout << data0(static_cast<int>(kEz), i, j);
  }
  return (time / repeats);
}

void RunAllModelsCalculation() {
  data1.resize(components, length_x + 2, length_y + 2);
  data0.resize(components, length_x + 2, length_y + 2);
  data1(kCeze, i, j) = one;
  data0(kCeze, i, j) = one;
  data1(kCezh, i, j) = coefficient;
  data0(kCezh, i, j) = coefficient;
  data1(kChxe, i, j) = coefficient;
  data0(kChxe, i, j) = coefficient;
  data1(kChxh, i, j) = coefficient;
  data0(kChxh, i, j) = coefficient;
  data1(kChye, i, j) = coefficient;
  data0(kChye, i, j) = coefficient;
  data1(kChyh, i, j) = coefficient;
  data0(kChyh, i, j) = coefficient;
  printf("\n>> __Processing task: 2D %dx%dx%d(steps)\n",
	 length_x, length_y, steps);
  UsingRange(length_x, length_y, steps);
  printf("\n>> __Done processing task: 2D %dx%dx%d(steps)\n",
	 length_x, length_y, steps);
}
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int total_processes_number, process_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  RunAllModelsCalculation();
  MPI_Finalize();
  return 0;
}

