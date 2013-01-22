/// @file   2D_test_std_0512_1600.c
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

#define SIZE 512
#define STEPS 1600
const int size = SIZE;
const int steps = STEPS;

const int depth = 2;
const int length_x = SIZE;
const int length_y = SIZE;

int sum_flag, time_flag, size_check_flag;
const double initial = 0.0, coefficient = 0.1;
enum DataComponents {kHx = 0, kChxh, kChxe,
                     kEz, kCeze, kCezh,
                     kHy, kChyh, kChye,
                     kSrcEz };
const int number_of_components = kSrcEz + 1;
const int d2 = SIZE + 2;
const int d1 = (SIZE + 2)*(SIZE + 2);
double *__restrict data0;
double *__restrict data1;
void UpdateFields2DStdC(int x1, int x2, int y1, int y2);
void InitializeFields();
void SetCoefficients();
void CycleArrays();
double EzControlSum();
double UsingStdC();

void initFDTD() {
  data0 = malloc(sizeof(double)*number_of_components*d1);
  data1 = malloc(sizeof(double)*number_of_components*d1);
  printf(">> __Now processing task 2D %dx%dx%d(steps)\n", length_x,
	 length_y, steps);
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
  printf(">> __Done processing task 2D %dx%dx%d(steps)\n", length_x, length_y,
       steps);
}
void UpdateFields2DStdC(int x1, int x2, int y1, int y2) {

  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
	data1[d1*kHx + d2*i3 + i4]
            = data0[d1*kChxh + d2*i3 + i4] *
            data0[d1*kHx + d2*i3 + i4]
            - data0[d1*kChxe + d2*i3 + i4] *
             (data0[d1*kEz + d2*i3 + i4 + 1]
              - data0[d1*kEz + d2*i3 + i4]);

	data1[d1*kHy + d2*i3 + i4]
            = data0[d1*kChyh + d2*i3 + i4] *
            data0[d1*kHy + d2*i3 + i4]
            + data0[d1*kChye + d2*i3 + i4] *
            (data0[d1*kEz + d2*(i3+1)  + i4]
              - data0[d1*kEz + d2*i3 + i4]);

	data1[d1*kEz + d2*i3 + i4]
            = data0[d1*kCeze + d2*i3 + i4] *
            data0[d1*kEz + d2*i3 + i4]
            + data0[d1*kCezh + d2*i3 + i4] *
	  ((data0[d1*kHy + d2*i3 + i4]
            - data0[d1*kHy + d2*(i3-1) + i4])
	   - (data0[d1*kHx + d2*i3 + i4]
              - data0[d1*kHx + d2*i3 + (i4-1)]))
            + data0[d1*kSrcEz + d2*i3 + i4];
    }
  }
}
void InitializeFields() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
        data0[d1*kHx + d2*i3 + i4] = initial;
        data0[d1*kHy + d2*i3 + i4] = initial;
        data0[d1*kEz + d2*i3 + i4] = initial;
        data0[d1*kSrcEz + d2*i3 + i4] = initial;
        data1[d1*kHx + d2*i3 + i4] = initial;
        data1[d1*kHy + d2*i3 + i4] = initial;
        data1[d1*kEz + d2*i3 + i4] = initial;
        data1[d1*kSrcEz + d2*i3 + i4] = initial;
    }
  }
}

void SetCoefficients() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
        data0[d1*kChxh + d2*i3 + i4] = coefficient;
        data0[d1*kChxe + d2*i3 + i4] = coefficient;
        data0[d1*kChyh + d2*i3 + i4] = coefficient;
        data0[d1*kChye + d2*i3 + i4] = coefficient;
        data0[d1*kCeze + d2*i3 + i4] = 1;
        data0[d1*kCezh + d2*i3 + i4] = coefficient;

        data1[d1*kChxh + d2*i3 + i4] = coefficient;
        data1[d1*kChxe + d2*i3 + i4] = coefficient;
        data1[d1*kChyh + d2*i3 + i4] = coefficient;
        data1[d1*kChye + d2*i3 + i4] = coefficient;
        data1[d1*kCeze + d2*i3 + i4] = 1;
        data1[d1*kCezh + d2*i3 + i4] = coefficient;
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
	sum += ABS(data0[d1*kEz + d2*i3 + i4]);
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
    src = exp(-1.0 * (double) (t * t) / (steps * steps));
    for (int i3 = 0; i3 <= length_x + 1; ++i3) {
      for (int i4 = 0; i4 <= length_y + 1; ++i4) {
          data0[d1*kSrcEz + d2*i3 + i4] = src;
      }
    }
    UpdateFields2DStdC(1, length_x, 1, length_y);
    CycleArrays();
  }
  clock_t finish_time = clock();
  double total_time = ((double)(finish_time - start_time))/CLOCKS_PER_SEC;
  if (time_flag == 1) printf(">> StdC 2D took %f\n", total_time);
  if (sum_flag == 1)
    printf(">> StdC 2D Ez control sum is %e\n",  EzControlSum());
  return total_time;
}

void RunCalculation() {
  initFDTD();
  UsingStdC();
  freeFDTD();
}

int main(int argc, char *argv[]) {
  // Flags *********************************************************
  sum_flag = 0;
  time_flag = 1;
  size_check_flag = 0;
  // ***************************************************************
  //printf("Mod2\n");
  RunCalculation();
  return 0;
}

