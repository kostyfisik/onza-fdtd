/// @file   3D_test_range_blocks.cc
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
#define BLOCK_SIZE sz_bl

// Setting model parameters
const int size = SIZE;
const int steps = STEPS;
const int length_x = SIZE;
const int length_y = SIZE;
const int length_z = SIZE;
const int block_size = BLOCK_SIZE;
const int repeats = 2;
const int sum_flag = 0, time_flag = 1, block_flag = bl_fl, print_flag = 0;
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

double RunBlocksRange(int length_x, int length_y, int length_z, int steps,
                            int block_size) {
  if (block_size >= size) {
    printf("\n>> RunBlocksRange:\n");
    printf("\n>> Error! Block size exceeds model size!\n");
    return 0.0;
  }
  else {
    int max = (length_x <= length_y) ? length_y:length_x;
    max = (length_z <= max) ? max:length_z;
    const int lim = ceil(max / block_size);
    const int sm_x_ = length_x / lim, sm_y_ = length_y / lim,
      sm_z_ = length_z / lim;
    double time = 0, start, end;
    const blitz::Range x_block_first(1, sm_x_), y_block_first(1, sm_y_),
      z_block_first(1, sm_z_);
    blitz::Array<blitz::Range, 1> Block_ranges_x(lim), Block_ranges_y(lim),
      Block_ranges_z(lim);
    for (int v = 1; v <= lim; ++v) {
      Block_ranges_x(v - 1) = x_block_first + (sm_x_ * (v - 1));
      Block_ranges_y(v - 1) = y_block_first + (sm_y_ * (v - 1));
      Block_ranges_z(v - 1) = z_block_first + (sm_z_ * (v - 1));
    }
    for (int l = 0; l < repeats; ++l) {
      data0(kHx, i, j, k) = initial;
      data0(kHy, i, j, k) = initial;
      data0(kHz, i, j, k) = initial;
      data0(kEx, i, j, k) = initial;
      data0(kEy, i, j, k) = initial;
      data0(kEz, i, j, k) = initial;
      start = MPI_Wtime();
      for (int t = 0; t < steps; ++t) {
        const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
        data0(static_cast<int>(kSrcEz), i, j, k) = src;
        for (int v1 = 1; v1 <= lim; ++v1) {
          for (int v2 = 1; v2 <= lim; ++v2) {
            for (int v3 = 1; v3 <= lim; ++v3) {
              UpdateFields3DRange(Block_ranges_x(v1 - 1),
                                  Block_ranges_y(v2 - 1),
                                  Block_ranges_z(v3 - 1));
            }
          }
        }
        cycleArrays(data0, data1);
      }
      end = MPI_Wtime();
      time += end - start;
    }
  if (time_flag == 1) {
    printf("\n>> Range 3D (%d blocks) took %f \n", lim,
           static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    printf("\n>> Range 3D (%d blocks) control sum for Ez is ", lim);
    std::cout << sum(abs(data0(static_cast<int>(kEz), i, j, k))) << "\n";
  }
  if (print_flag == 1) {
    std::cout << data0(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
  }
}

void RunAllModelsCalculation() {
  if (block_flag == 1) {
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
    RunBlocksRange(length_x, length_y, length_z, steps, block_size);
  printf("\n>> __Done processing task: 3D %dx%dx%dx%d(steps)\n",
	 length_x, length_y, length_z, steps);
  }
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

