/// @file   2D_test_stencil.cc
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
#define CONTROL_SUM cs

// Setting model parameters
const int size = SIZE;
const long int steps = STEPS;
const int length_x = SIZE;
const int length_y = SIZE;
const int repeats = 1;
const int sum_flag = CONTROL_SUM, time_flag = 1, print_flag = 0;
const double initial = 0.0, coefficient = 0.1, one = 1.0;
// End Setting model parameters

blitz::Array <double, 3> data0, data1;
const blitz::Range i = blitz::Range::all(), j = blitz::Range::all();

enum DataComponents {kEz = 0, kCeze, kCezh,
                     kHx, kChxh, kChxe,
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
    blitz::Array<double, 2> tmpkHx = (data1(kHx, i, j1));
    const blitz::Range i1(-1 + i(0), i(i.length() - 1));
    blitz::Array<double, 2> tmpkHy = (data1(kHy, i1, j));
    applyStencil(update_field_Ez(), tmpkEz_next, tmpkEz, tmpkCeze,
                 tmpkCezh, tmpkHy, tmpkHx, tmpkSrcEz);
    // cycleArrays(data1(kEz, i, j), tmpkEz_next);
       data1(kEz, i, j) = tmpkEz_next;
  }
}
double UsingStencils(int length_x, int length_y, long int steps) {
  double time = 0, start, end;
  const blitz::Range first(1, length_x), second(1, length_y);
  for (int o = 0; o < repeats; ++o) {
    data0(kEz, i, j) = initial;
    data0(kHx, i, j) = initial;
    data0(kHy, i, j) = initial;
    data1(kEz, i, j) = initial;
    data1(kHx, i, j) = initial;
    data1(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (long int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double> (t) * static_cast<double> (t)
                               / (static_cast<double> (steps) * static_cast<double> (steps)));
      data0(static_cast<int>(kSrcEz), i, j) = src;
      UpdateFields2DStencil(first, second);
      cycleArrays(data0, data1);
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf(">> 2D Stencil %dx%dx%ld(steps) took %f \n",
           length_x, length_y, steps, static_cast<float>(time / repeats));
  }
  if (cs == 1) {
  std::cout << ">> 3D Stencil control sum for Ez is ";
  printf("%23.16e\n", sum(data0(static_cast<int>(kEz), i, j)));
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
  printf("\n>> 2D Stencil %04d\n",
	 length_x);
  fflush(stdout);
  UsingStencils(length_x, length_y, steps);
  // printf("\n>> __Done processing task: 2D %dx%dx%d(steps)\n",
  //        length_x, length_y, steps);
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

