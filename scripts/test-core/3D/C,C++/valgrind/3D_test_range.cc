/// @file   3D_test_range.cc
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
const int length_z = SIZE;
const int repeats = 2;
const int sum_flag = 0, time_flag = 1, print_flag = 0;
const double initial = 0.0, coefficient = 0.1, one = 1.0;
// End Setting model parameters

blitz::Array <double, 4> data0, data1;
const blitz::Range i = blitz::Range::all(), j = blitz::Range::all(),
  k = blitz::Range::all();

enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze,
                     kSrcEz };
const int components = kSrcEz + 1;

void UpdateFields3DRange(blitz::Range i, blitz::Range j, blitz::Range k) {
  data1(kHx, i, j, k) = data0(kChxh, i, j, k) *
          data0(kHx, i, j, k) + data0(kChxe, i, j, k) *
          ((data0(kEy, i, j, k + 1) - data0(kEy, i, j, k))
           - (data0(kEz, i, j + 1, k) - data0(kEz, i, j, k)));
 
  data1(kHy, i, j, k) = data0(kChyh, i, j, k) *
          data0(kHy, i, j, k) + data0(kChye, i, j, k) *
          ((data0(kEz, i + 1, j, k) - data0(kEz, i, j, k))
           - (data0(kEx, i, j, k + 1) - data0(kEx, i, j, k)));
  
  data1(kHz, i, j, k) = data0(kChzh, i, j, k) *
          data0(kHz, i, j, k) + data0(kChze, i, j, k) *
          ((data0(kEx, i, j + 1, k) - data0(kEx, i, j, k))
           - (data0(kEy, i + 1, j, k) - data0(kEy, i, j, k)));
  
  data1(kEx, i, j, k) = data0(kCexe, i, j, k) *
          data0(kEx, i, j, k) + data0(kCexh, i, j, k) *
          ((data0(kHz, i, j, k) - data0(kHz, i, j - 1, k))
           - (data0(kHy, i, j, k) - data0(kHy, i, j, k - 1)));
  
  data1(kEy, i, j, k) = data0(kCeye, i, j, k) *
          data0(kEy, i, j, k) + data0(kCeyh, i, j, k) *
          ((data0(kHx, i, j, k) - data0(kHx, i, j, k - 1))
           - (data0(kHz, i, j, k) - data0(kHy, i - 1, j, k)));

  data1(kEz, i, j, k) = data0(kCeze, i, j, k) *
          data0(kEz, i, j, k) + data0(kCezh, i, j, k) *
          ((data0(kHy, i, j,  k) - data0(kHy, i - 1, j, k))
           - (data0(kHx, i, j, k) - data0(kHx, i, j - 1, k)))
          + data0(kSrcEz, i, j, k);
}

double UsingRange(int length_x, int length_y, int length_z, int steps) {
  const blitz::Range inner_i(1, length_x);
  const blitz::Range inner_j(1, length_y);
  const blitz::Range inner_k(1, length_z);
  double time = 0, start, end;
  for (int l = 0; l < repeats; ++l) {
    data0(kHx, i, j, k) = initial;
    data0(kHy, i, j, k) = initial;
    data0(kHz, i, j, k) = initial;
    data0(kEx, i, j, k) = initial;
    data0(kEy, i, j, k) = initial;
    data0(kEz, i, j, k) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src = exp(-static_cast<double>(t * t) / (steps * steps));
      data0(static_cast<int>(kSrcEz), i, j, k) = src;
      UpdateFields3DRange(inner_i, inner_j, inner_k);
      cycleArrays(data0, data1);
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Range 3D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Range 3D control sum for Ez is "
              << sum(abs(data0(static_cast<int>(kEz), i, j, k))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data0(static_cast<int>(kHx), i, j, k);
    std::cout << data0(static_cast<int>(kHy), i, j, k);
    std::cout << data0(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
}

void RunAllModelsCalculation() {
  data1.resize(components, length_x + 2, length_y + 2, length_z + 2);
  data0.resize(components, length_x + 2, length_y + 2, length_z + 2);
  data1(kCeze, i, j, k) = one;
  data0(kCeze, i, j, k) = one;
  data1(kCexe, i, j, k) = coefficient;
  data0(kCexe, i, j, k) = coefficient;
  data1(kCexh, i, j, k) = coefficient;
  data0(kCexh, i, j, k) = coefficient;
  data1(kCeye, i, j, k) = coefficient;
  data0(kCeye, i, j, k) = coefficient;
  data1(kCeyh, i, j, k) = coefficient;
  data0(kCeyh, i, j, k) = coefficient;
  data1(kCezh, i, j, k) = coefficient;
  data0(kCezh, i, j, k) = coefficient;
  data1(kChxe, i, j, k) = coefficient;
  data0(kChxe, i, j, k) = coefficient;
  data1(kChxh, i, j, k) = coefficient;
  data0(kChxh, i, j, k) = coefficient;
  data1(kChye, i, j, k) = coefficient;
  data0(kChye, i, j, k) = coefficient;
  data1(kChyh, i, j, k) = coefficient;
  data0(kChyh, i, j, k) = coefficient;
  data1(kChzh, i, j, k) = coefficient;
  data0(kChzh, i, j, k) = coefficient;
  data1(kChze, i, j, k) = coefficient;
  data0(kChze, i, j, k) = coefficient;
  printf("\n>> __Processing task: 3D %dx%dx%dx%d(steps)\n",
	 length_x, length_y, length_z, steps);
  UsingRange(length_x, length_y, length_z, steps);
  printf("\n>> __Done processing task: 3D %dx%dx%dx%d(steps)\n",
	 length_x, length_y, length_z, steps);
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

