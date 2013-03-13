/// @file   2D_test_stencil_blocks.cc
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
const int depth = 2;
const int size = SIZE;
const int steps = STEPS;
const int length_x = SIZE;
const int length_y = SIZE;
const int block_size = BLOCK_SIZE;
const int repeats = 2;
const int sum_flag = 0, time_flag = 1, block_flag = bl_fl,
  print_flag = 0;
const double initial = 0.0, coefficient = 0.1, one = 1.0;
// End setting model parameters

blitz::Array <double, 3> data0, data1;
const blitz::Range i = blitz::Range::all(), j = blitz::Range::all();
enum DataComponents {kHx = 0, kChxh, kChxe,
                     kEz, kCeze, kCezh,
                     kHy, kChyh, kChye,
                     kSrcEz};
const int components = kSrcEz + 1;

BZ_DECLARE_STENCIL5(update_field_Hx, kHx_next, kHx, kChxh, kChxe, kEz)
kHx_next = kChxh * kHx - kChxe * (kEz(0, 1) - kEz);
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, 0), blitz::shape(0, 1))

BZ_DECLARE_STENCIL5(update_field_Hy, kHy_next, kHy, kChyh, kChye, kEz)
kHy_next = kChyh * kHy + kChye * (kEz(1, 0) - kEz);
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, 0), blitz::shape(1, 0))

BZ_DECLARE_STENCIL7(update_field_Ez, kEz_next, kEz, kCeze, kCezh,
                           kHy, kHx, kSrcEz)
kEz_next = kCeze * kEz + kCezh * ((kHy(1, 0) - kHy) - (kHx(0, 1) - kHx));
kEz_next = kEz_next + kSrcEz;
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(1, 1), blitz::shape(0, 0))

void UpdateFields2DStencil(blitz::Range i, blitz::Range j) {
  {
    const blitz::Range j1(j(0), 1 + j(j.length() - 1)); 
    blitz::Array<double, 2> tmpkHx_next = (data1(kHx, i, j1));
    blitz::Array<double, 2> tmpkHx = (data0(kHx, i, j1));
    blitz::Array<double, 2> tmpkChxh = (data0(kChxh, i, j1));
    blitz::Array<double, 2> tmpkChxe = (data0(kChxe, i, j1));
    blitz::Array<double, 2> tmpkEz = (data0(kEz, i, j1));
    applyStencil(update_field_Hx(), tmpkHx_next, tmpkHx, tmpkChxh,
                 tmpkChxe, tmpkEz);
    //    cycleArrays(data1(kHx, i, j1), tmpkHx_next);
      data1(kHx, i, j1) = tmpkHx_next;
  }
  {
    const blitz::Range i1(i(0), 1 + i(i.length() - 1));
    blitz::Array<double, 2> tmpkHy_next = (data1(kHy, i1, j));
    blitz::Array<double, 2> tmpkHy = (data0(kHy, i1, j));
    blitz::Array<double, 2> tmpkChyh = (data0(kChyh, i1, j));
    blitz::Array<double, 2> tmpkChye = (data0(kChye, i1, j));
    blitz::Array<double, 2> tmpkEz = (data0(kEz, i1, j));
    applyStencil(update_field_Hy(), tmpkHy_next, tmpkHy, tmpkChyh,
                 tmpkChye, tmpkEz);
    // cycleArrays(data1(kHy, i1, j), tmpkHy_next);
       data1(kHy, i1, j) = tmpkHy_next;
  }
  {
    blitz::Array<double, 2> tmpkEz_next = (data1(kEz, i, j));
    blitz::Array<double, 2> tmpkEz = (data0(kEz, i, j));
    blitz::Array<double, 2> tmpkCeze = (data0(kCeze, i, j));
    blitz::Array<double, 2> tmpkCezh = (data0(kCezh, i, j));
    blitz::Array<double, 2> tmpkSrcEz = (data0(kSrcEz, i, j));
    const blitz::Range j1(-1 + j(0), j(j.length() - 1));
    blitz::Array<double, 2> tmpkHx = (data0(kHx, i, j1));
    const blitz::Range i1(-1 + i(0), i(i.length() - 1));
    blitz::Array<double, 2> tmpkHy = (data0(kHy, i1, j));
    applyStencil(update_field_Ez(), tmpkEz_next, tmpkEz, tmpkCeze,
                 tmpkCezh, tmpkHy, tmpkHx, tmpkSrcEz);
    // cycleArrays(data1(kEz, i, j), tmpkEz_next);
       data1(kEz, i, j) = tmpkEz_next;
  }
}

double RunBlocksStencil(int length_x, int length_y, int steps,
                            int block_size) {
  if (block_size >= size) {
    printf("\n>> RunBlocksStencil:\n");
    printf("\n>> Error! Block size exceeds model size!\n");
    return 0.0;
  }
  else {
    int blocks_number_y = ceil(length_y / block_size);
    int blocks_number_x = ceil(length_x / block_size);
    int lim, sm_x, sm_y;
    if (blocks_number_x > blocks_number_y) {
      lim = blocks_number_x;
      sm_x = block_size;
      sm_y = length_y / lim;
    } else {
      lim = blocks_number_y;
      sm_y = block_size;
      sm_x = length_x / lim;
    }
    const int sm_x_ = sm_x;
    const int sm_y_ = sm_y;
    const int lim_ = lim;
    double time = 0, start, end;
    const blitz::Range x_block_first(1, sm_x_), y_block_first(1, sm_y_);
    blitz::Array<blitz::Range, 1> Block_ranges_x(lim_), Block_ranges_y(lim_);
    for (int v = 1; v <= lim_; ++v) {
      Block_ranges_x(v - 1) = x_block_first + (sm_x_ * (v - 1));
      Block_ranges_y(v - 1) = y_block_first + (sm_y_ * (v - 1));
    }
    for (int l = 0; l < repeats; ++l) {
      data0(kEz, i, j) = initial;
      data0(kHx, i, j) = initial;
      data0(kHy, i, j) = initial;
      start = MPI_Wtime();
      for (int t = 0; t < steps; ++t) {
        const double src = exp(-static_cast<double>(t * t) / (steps * steps));
        data0(static_cast<int>(kSrcEz), i, j) = src;
        for (int v1 = 1; v1 <= lim_; ++v1) {
          for (int v2 = 1; v2 <= lim_; ++v2) {
            UpdateFields2DStencil(Block_ranges_x(v1 - 1), Block_ranges_y(v2 - 1));
          }
        }
        cycleArrays(data0, data1);
      }
      end = MPI_Wtime();
      time += end - start;
    }
    if (time_flag == 1) {
      printf("\n>> Stencil 2D(%d blocks) took %f \n", lim_,
             static_cast<float>(time / repeats));
    }
    if (sum_flag == 1) {
      printf("\n>> Stencil 2D (%d blocks) control sum for Ez is ", lim_);
      std::cout << sum(abs(data0(static_cast<int>(kEz), i, j))) << std::endl;
    }
    if (print_flag == 1) {
      std::cout << data0(static_cast<int>(kHx), i, j);
      std::cout << data0(static_cast<int>(kHy), i, j);
      std::cout << data0(static_cast<int>(kEz), i, j);
    }
    return (time / repeats);
  }
}

void RunAllModelsCalculation() {
  if (block_flag == 1) {
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
    RunBlocksStencil(length_x, length_y, steps, block_size);
    printf("\n>> __Done processing task: 2D %dx%dx%d(steps)\n",
           length_x, length_y, steps);
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

