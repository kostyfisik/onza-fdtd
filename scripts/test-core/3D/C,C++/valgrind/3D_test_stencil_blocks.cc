/// @file   3D_test_stencil_blocks.cc
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

BZ_DECLARE_STENCIL6(update_field_Hx, kHx_next, kChxh, kHx, kChxe, kEy, kEz)
kHx_next = kChxh * kHx + kChxe * ((kEy(0, 0, 1) - kEy) - (kEz(0, 1, 0) - kEz));
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, 0, 0), blitz::shape(0, 1, 1))

BZ_DECLARE_STENCIL6(update_field_Hy, kHy_next, kChyh, kHy, kChye, kEz, kEx)
kHy_next = kChyh * kHy + kChye * ((kEz(1, 0, 0) - kEz) - (kEx(0, 0, 1) - kEx));
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, 0, 0), blitz::shape(1, 0, 1))

BZ_DECLARE_STENCIL6(update_field_Hz, kHz_next, kChzh, kHz, kChze, kEx, kEy)
kHz_next = kChzh * kHz + kChze * ((kEx(0, 1, 0) - kEx) - (kEy(1, 0, 0) - kEy));
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, 0, 0), blitz::shape(1, 1, 0))

BZ_DECLARE_STENCIL6(update_field_Ex, kEx_next, kCexe, kEx, kCexh, kHz, kHy)
kEx_next = kCexe * kEx + kCexh *
  ((kHz - kHz(0, -1, 0)) - (kHy - kHy(0, 0, -1)));
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(0, -1, -1), blitz::shape(0, 0, 0))

BZ_DECLARE_STENCIL6(update_field_Ey, kEy_next, kCeye, kEy, kCeyh, kHx, kHz)
kEy_next = kCeye * kEy + kCeyh *
  ((kHx - kHx(0, 0, -1)) - (kHz - kHz(-1, 0, 0)));
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(-1, 0, -1), blitz::shape(0, 0, 0))

BZ_DECLARE_STENCIL7(update_field_Ez, kEz_next, kCeze, kEz, kCezh,
                    kHy, kHx, kSrcEz)
kEz_next = kCeze * kEz + kCezh *
  ((kHy - kHy(-1, 0, 0)) - (kHx - kHx(0, -1, 0)));
kEz_next = kEz_next + kSrcEz;
BZ_END_STENCIL_WITH_SHAPE(blitz::shape(-1, -1, 0), blitz::shape(0, 0, 0))

void UpdateFields3DStencil(blitz::Range i, blitz::Range j, blitz::Range k) {
 
  {
    const blitz::Range j1(j(0), 1 + j(j.length() - 1)),
      k1(k(0), 1 + k(k.length() - 1));
    blitz::Array<double, 3> tmpkHx_next = (data1(kHx, i, j1, k1));
    blitz::Array<double, 3> tmpkChxh = (data0(kChxh, i, j1, k1));
    blitz::Array<double, 3> tmpkHx = (data0(kHx, i, j1, k1));
    blitz::Array<double, 3> tmpkChxe = (data0(kChxe, i, j1, k1));
    blitz::Array<double, 3> tmpkEy = (data0(kEy, i, j1, k1));
    blitz::Array<double, 3> tmpkEz = (data0(kEz, i, j1, k1));
    applyStencil(update_field_Hx(), tmpkHx_next, tmpkChxh, tmpkHx,
                 tmpkChxe, tmpkEy, tmpkEz);
    data1(kHx, i, j1, k1) = tmpkHx_next;
  }
  {
    const blitz::Range i1(i(0), 1 + i(i.length() - 1)),
      k1(k(0), 1 + k(k.length() - 1));
    blitz::Array<double, 3> tmpkHy_next = (data1(kHy, i1, j, k1));
    blitz::Array<double, 3> tmpkChyh = (data0(kChyh, i1, j, k1));
    blitz::Array<double, 3> tmpkHy = (data0(kHy, i1, j, k1));
    blitz::Array<double, 3> tmpkChye = (data0(kChye, i1, j, k1));
    blitz::Array<double, 3> tmpkEz = (data0(kEz, i1, j, k1));
    blitz::Array<double, 3> tmpkEx = (data0(kEx, i1, j, k1));
    applyStencil(update_field_Hy(), tmpkHy_next, tmpkChyh, tmpkHy,
                 tmpkChye, tmpkEz, tmpkEx);
    data1(kHy, i1, j, k1) = tmpkHy_next;
  }
  {
  const blitz::Range i1(i(0), 1 + i(i.length() - 1)),
      j1(j(0), 1 + j(j.length() - 1));
    blitz::Array<double, 3> tmpkHz_next = (data1(kHz, i1, j1, k));
    blitz::Array<double, 3> tmpkChzh = (data0(kChzh, i1, j1, k));
    blitz::Array<double, 3> tmpkHz = (data0(kHz, i1, j1, k));
    blitz::Array<double, 3> tmpkChze = (data0(kChze, i1, j1, k));
    blitz::Array<double, 3> tmpkEx = (data0(kEx, i1, j1, k));
    blitz::Array<double, 3> tmpkEy = (data0(kEy, i1, j1, k));
    applyStencil(update_field_Hz(), tmpkHz_next, tmpkChzh, tmpkHz,
                 tmpkChze, tmpkEx, tmpkEy);
    data1(kHz, i1, j1, k) = tmpkHz_next;
  }
  {
    const blitz::Range j1(-1 + j(0), j(j.length() - 1)),
      k1(-1 + k(0), k(k.length() - 1));
    blitz::Array<double, 3> tmpkEx_next = (data1(kEx, i, j1, k1));
    blitz::Array<double, 3> tmpkCexe = (data0(kCexe, i, j1, k1));
    blitz::Array<double, 3> tmpkEx = (data0(kEx, i, j1, k1));
    blitz::Array<double, 3> tmpkCexh = (data0(kCexh, i, j1, k1));
    blitz::Array<double, 3> tmpkHz = (data0(kHz, i, j1, k1));
    blitz::Array<double, 3> tmpkHy = (data0(kHy, i, j1, k1));
    applyStencil(update_field_Ex(), tmpkEx_next, tmpkCexe, tmpkEx,
                 tmpkCexh, tmpkHz, tmpkHy);
    data1(kEx, i, j1, k1) = tmpkEx_next;
  }
  {
    const blitz::Range i1(-1 + i(0), i(i.length() - 1)),
      k1(-1 + k(0), k(k.length() - 1));
    blitz::Array<double, 3> tmpkEy_next = (data1(kEy, i1, j, k1));
    blitz::Array<double, 3> tmpkCeye = (data0(kCeye, i1, j, k1));
    blitz::Array<double, 3> tmpkEy = (data0(kEy, i1, j, k1));
    blitz::Array<double, 3> tmpkCeyh = (data0(kCeyh, i1, j, k1)); 
    blitz::Array<double, 3> tmpkHx = (data0(kHx, i1, j, k1));
    blitz::Array<double, 3> tmpkHz = (data0(kHz, i1, j, k1));
    applyStencil(update_field_Ey(), tmpkEy_next, tmpkCeye, tmpkEy,
                 tmpkCeyh, tmpkHx, tmpkHz);
    data1(kEy, i1, j, k1) = tmpkEy_next;
  }
  {
    const blitz::Range i1(-1 + i(0), i(i.length() - 1)),
      j1(-1 + j(0), j(j.length() - 1));  
    blitz::Array<double, 3> tmpkEz_next = (data1(kEz, i1, j1, k));
    blitz::Array<double, 3> tmpkCeze = (data0(kCeze, i1, j1, k));
    blitz::Array<double, 3> tmpkEz = (data0(kEz, i1, j1, k));
    blitz::Array<double, 3> tmpkCezh = (data0(kCezh, i1, j1, k));
    blitz::Array<double, 3> tmpkSrcEz = (data0(kSrcEz, i1, j1, k));
    blitz::Array<double, 3> tmpkHy = (data0(kHy, i1, j1, k));
    blitz::Array<double, 3> tmpkHx = (data0(kHx, i1, j1, k));
    applyStencil(update_field_Ez(), tmpkEz_next, tmpkCeze, tmpkEz,
                 tmpkCezh, tmpkHy, tmpkHx, tmpkSrcEz);
    data1(kEz, i1, j1, k) = tmpkEz_next;
  }
}

double RunBlocksStencil(int length_x, int length_y, int length_z,
				   int steps, int block_size) {
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
              UpdateFields3DStencil(Block_ranges_x(v1 - 1),
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
      printf("\n>> Stencil 3D (%d blocks) took %f \n", lim,
             static_cast<float>(time / repeats));
    }
    if (sum_flag == 1) {
      printf("\n>> Stencil 3D (%d blocks) control sum for Ez is ", lim);
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
    RunBlocksStencil(length_x, length_y, length_z, steps, block_size);
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
