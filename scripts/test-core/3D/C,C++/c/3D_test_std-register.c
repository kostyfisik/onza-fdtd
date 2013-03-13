/// @file   3D_test_c_register.cc
/// @author Markovich Dmitry <dmmrkovich at gmail (.) com>
/// @copyright 2012 Markovich Dmitry
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2013 Ladutenko Konstantin
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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <xmmintrin.h>

#define SIZE sz
#define STEPS st

const int size = SIZE;
const int steps = STEPS;

const int depth = 2;
const int length_x = SIZE;
const int length_y = SIZE;
const int length_z = SIZE;

int sum_flag, time_flag, size_check_flag;
const double initial = 0.0, coefficient = 0.1;
enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze, 
                     kSrcEz };
const int number_of_components = kSrcEz + 1;
const int s_y = SIZE + 2;
const int s_x = (SIZE + 2)*(SIZE + 2);
const int dc = (SIZE + 2)*(SIZE + 2)*(SIZE + 2);
double *__restrict data0;
double *__restrict data1;
void UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2, int t);
void InitializeFields();
void SetCoefficients();
void CycleArrays();
double EzControlSum();
double UsingStdC();

void initFDTD() {
  data0 = malloc(sizeof(double)*number_of_components*dc);
  data1 = malloc(sizeof(double)*number_of_components*dc);
  printf("\n>> __Now processing task 3D %dx%dx%dx%d(steps)\n", length_x,
	 length_y, length_z, steps);
}
double ABS(double number) {
  if (number >= 0) {
    return number;
  } else {
    return -1 * number;
  }
}
void  freeFDTD() {
  free(data0);
  free(data1);
  printf("\n>> _ _Done processing task 3D %dx%dx%dx%d(steps)\n", length_x, length_y,
       length_z, steps);
}
void UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2, int t) {
  double t1kHx, t0kHx, t1kHy, t0kHy, t1kHz, t0kHz,
      t1kEx, t0kEx, t1kEy, t0kEy, t1kEz, t0kEz,
      t0kChxh, t0kChxe, t0kChyh, t0kChye, t0kChzh, t0kChze,
      t0kCexh, t0kCexe, t0kCeyh, t0kCeye, t0kCezh, t0kCeze,
      t0kEx1py, t0kEx1pz, t0kEy1px, t0kEy1pz, t0kEz1px, t0kEz1py,
      t0kHx1my, t0kHx1mz, t0kHy1mx, t0kHy1mz, t0kHz1mx, t0kHz1my,
      t0kSrcEz;
  double src = exp(-1.0 * (double) (t * t) / (steps * steps));
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; i5+=1) {       
        const int cell = s_x*i3 + s_y*i4 + i5;
        data0[dc*kSrcEz + cell] = src;
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      t0kEy = data0[dc*kEy + s_x*i3 + s_y*i4 + z1];
      t0kEx = data0[dc*kEx + s_x*i3 + s_y*i4 + z1];
      t0kHx1mz = data0[dc*kHx + s_x*i3 + s_y*i4 + z1 - 1];
      t0kHy1mz = data0[dc*kHy + s_x*i3 + s_y*i4 + z1 - 1];
      for (int i5 = z1; i5 <= z2; i5+=1) {       
        const int cell = s_x*i3 + s_y*i4 + i5;
        t0kHx = data0[dc*kHx + cell];
        t0kHy = data0[dc*kHy + cell];
        t0kHz = data0[dc*kHz + cell];
        t0kEy1pz = data0[dc*kEy + cell+1];
        t0kEx1pz = data0[dc*kEx + cell+1];
        t0kEz = data0[dc*kEz + cell];
        data1[dc*kHx + cell]
            = data0[dc*kChxh + cell] * t0kHx
            + data0[dc*kChxe + cell] *
            ( (t0kEy1pz - t0kEy) - (data0[dc*kEz + cell + s_y] - t0kEz));
        data1[dc*kHy + cell]
            = data0[dc*kChyh + cell] * t0kHy
            + data0[dc*kChye + cell] *
            ((data0[dc*kEz + cell + s_x]- t0kEz) - (t0kEx1pz - t0kEx));
	data1[dc*kHz + cell]
            = data0[dc*kChzh + cell] * t0kHz
            + data0[dc*kChze + cell] *
            ((data0[dc*kEx + cell + s_y] - t0kEx)
             - (data0[dc*kEy + cell + s_x] - t0kEy));
	data1[dc*kEx + cell]
            = data0[dc*kCexe + cell] * t0kEx
            + data0[dc*kCexh + cell] *
            ((t0kHz - data0[dc*kHz + cell - s_y]) - (t0kHy - t0kHy1mz));
        data1[dc*kEy + cell]
            = data0[dc*kCeye + cell] * t0kEy 
            + data0[dc*kCeyh + cell] *
            ((t0kHx - t0kHx1mz) - (t0kHz - data0[dc*kHz + cell - s_x]));
	data1[dc*kEz + cell]
            = data0[dc*kCeze + cell] * t0kEz
            + data0[dc*kCezh + cell] *
            ((t0kHy - data0[dc*kHy + cell - s_x])
             - (t0kHx - data0[dc*kHx + cell - s_y]))
            + data0[dc*kSrcEz + cell];
        t0kEy = t0kEy1pz;
        t0kEx = t0kEx1pz;
        t0kHx1mz = t0kHx;
        t0kHy1mz = t0kHy;
      }            
    }
  }
}
void InitializeFields() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z + 1; ++i5) {
        data0[dc*kHx + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kHy + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kHz + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kEx + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kEy + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kEz + s_x*i3 + s_y*i4 + i5] = initial;
        data0[dc*kSrcEz + s_x*i3 + s_y*i4 + i5] = initial;

        data1[dc*kHx + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kHy + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kHz + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kEx + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kEy + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kEz + s_x*i3 + s_y*i4 + i5] = initial;
        data1[dc*kSrcEz + s_x*i3 + s_y*i4 + i5] = initial;
      }
    }
  }
}

void SetCoefficients() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z + 1; ++i5) {
        data0[dc*kChxh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kChxe + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kChyh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kChye + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kChzh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kChze + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kCexe + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kCexh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kCeye + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kCeyh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data0[dc*kCeze + s_x*i3 + s_y*i4 + i5] = 1;
        data0[dc*kCezh + s_x*i3 + s_y*i4 + i5] = coefficient;

        data1[dc*kChxh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kChxe + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kChyh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kChye + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kChzh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kChze + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kCexe + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kCexh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kCeye + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kCeyh + s_x*i3 + s_y*i4 + i5] = coefficient;
        data1[dc*kCeze + s_x*i3 + s_y*i4 + i5] = 1;
        data1[dc*kCezh + s_x*i3 + s_y*i4 + i5] = coefficient;
      }
    }
  } 
}

void CycleArrays() {
  double *tmp;
  tmp = data0;
  data0 = data1;
  data1 = tmp;
}

double EzControlSum() {
  double sum = 0;
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z + 1; ++i5) {
	sum += ABS(data0[dc*kEz + s_x*i3 + s_y*i4 + i5]);
      }
    }
  }
  return sum;
}

double UsingStdC() {
  SetCoefficients();
  InitializeFields();
  clock_t start_time = clock();
  double src;
  for (int t = 0; t < steps; ++t) {
    UpdateFields3DStdC(1, length_x, 1, length_y, 1, length_z, t);
    CycleArrays();
  }
  clock_t finish_time = clock();

  double total_time = ((double)(finish_time - start_time))/CLOCKS_PER_SEC;
  if (time_flag == 1) printf(">> StdC 3D took %f\n", total_time);
  if (sum_flag == 1)
    printf("\n>> StdC 3D Ez control sum is %11.9e\n",  EzControlSum());
  return total_time;
}

void RunCalculation() {
  initFDTD();
  UsingStdC();
  freeFDTD();
}

int main(int argc, char *argv[]) {
  // Flags *********************************************************
  sum_flag = 1;
  time_flag = 1;
  size_check_flag = 0;
  // ***************************************************************
  //printf("Mod2\n");
  RunCalculation();
  return 0;
}

