/// @file   3D_test_c.cc
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
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
//#include <xmmintrin.h>

#define SIZE sz
#define STEPS st
#define CONTROL_SUM cs

const int size = SIZE;
const long int steps = STEPS;

const int depth = 2;
const int length_x = SIZE;
const int length_y = SIZE;
const int length_z = SIZE;
unsigned long long int num_of_updates = 0;
int sum_flag = CONTROL_SUM, time_flag = 1, size_check_flag = 1;
const double initial = 0.0, coefficient = 0.1;
enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze, 
                     kSrcEz };
const int number_of_components = kSrcEz + 1;
const int d3 = SIZE + 2;
const int d2 = (SIZE + 2)*(SIZE + 2);
const int d1 = (SIZE + 2)*(SIZE + 2)*(SIZE + 2);
double *__restrict data0;
double *__restrict data1;
void UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2);
void InitializeFields();
void SetCoefficients();
void CycleArrays();
double EzControlSum();
double UsingStdC();

void initFDTD() {
  data0 = malloc(sizeof(double)*number_of_components*d1);
  data1 = malloc(sizeof(double)*number_of_components*d1);
  printf("\n>> 3D StdC %04d\n",
	 length_x);
  fflush(stdout);
 // printf("\n>> __Now processing task 3D %dx%dx%dx%d(steps)\n", length_x,
//	 length_y, length_z, steps);
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
 // printf("\n>> __Done processing task 3D %dx%dx%dx%d(steps)\n", length_x, length_y,
 //      length_z, steps);
}
void UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2) {

  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      
      for (int i5 = z1; i5 <= z2; i5+=2) {
	data1[d1*kHx + d2*i3 + d3*i4 + i5]
            = data0[d1*kChxh + d2*i3 + d3*i4 + i5] *
            data0[d1*kHx + d2*i3 + d3*i4 + i5]
            + data0[d1*kChxe + d2*i3 + d3*i4 + i5] *
            ((data0[d1*kEy + d2*i3 + d3*i4 + i5+1]
              - data0[d1*kEy + d2*i3 + d3*i4 + i5])
             -(data0[d1*kEz + d2*i3 + d3*(i4+1)  + i5]
               - data0[d1*kEz + d2*i3 + d3*i4 + i5]));
	
	data1[d1*kHx + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kChxh + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kHx + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kChxe + d2*i3 + d3*i4 + i5+1] *
            ((data0[d1*kEy + d2*i3 + d3*i4 + i5+1+1]
              - data0[d1*kEy + d2*i3 + d3*i4 + i5+1])
             -(data0[d1*kEz + d2*i3 + d3*(i4+1)  + i5+1]
               - data0[d1*kEz + d2*i3 + d3*i4 + i5+1]));
	}
      
      for (int i5 = z1; i5 <= z2; i5+=2) {
	data1[d1*kHy + d2*i3 + d3*i4 + i5]
            = data0[d1*kChyh + d2*i3 + d3*i4 + i5] *
            data0[d1*kHy + d2*i3 + d3*i4 + i5]
            + data0[d1*kChye + d2*i3 + d3*i4 + i5] *
            ((data0[d1*kEz + d2*(i3+1)  + d3*i4 + i5]
              - data0[d1*kEz + d2*i3 + d3*i4 + i5])
             -(data0[d1*kEx + d2*i3 + d3*i4 + i5+1]
               - data0[d1*kEx + d2*i3 + d3*i4 + i5]));
        
        data1[d1*kHy + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kChyh + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kHy + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kChye + d2*i3 + d3*i4 + i5+1] *
            ((data0[d1*kEz + d2*(i3+1)  + d3*i4 + i5+1]
              - data0[d1*kEz + d2*i3 + d3*i4 + i5+1])
             -(data0[d1*kEx + d2*i3 + d3*i4 + i5+1+1]
               - data0[d1*kEx + d2*i3 + d3*i4 + i5+1]));
      }
      
      for (int i5 = z1; i5 <= z2; i5+=2) {
	data1[d1*kHz + d2*i3 + d3*i4 + i5]
            = data0[d1*kChzh + d2*i3 + d3*i4 + i5] *
            data0[d1*kHz + d2*i3 + d3*i4 + i5]
            + data0[d1*kChze + d2*i3 + d3*i4 + i5] *
            ((data0[d1*kEx + d2*i3 + d3*(i4+1) + i5]
              - data0[d1*kEx + d2*i3 + d3*i4 + i5])
             -(data0[d1*kEy + d2*(i3+1) + d3*i4 + i5]
             - data0[d1*kEy + d2*i3 + d3*i4 + i5]));

        data1[d1*kHz + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kChzh + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kHz + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kChze + d2*i3 + d3*i4 + i5+1] *
            ((data0[d1*kEx + d2*i3 + d3*(i4+1) + i5+1]
              - data0[d1*kEx + d2*i3 + d3*i4 + i5+1])
             -(data0[d1*kEy + d2*(i3+1) + d3*i4 + i5+1]
             - data0[d1*kEy + d2*i3 + d3*i4 + i5+1]));
      }
      
      for (int i5 = z1; i5 <= z2; i5+=2) {
	data1[d1*kEx + d2*i3 + d3*i4 + i5]
            = data0[d1*kCexe + d2*i3 + d3*i4 + i5] *
            data0[d1*kEx + d2*i3 + d3*i4 + i5]
            + data0[d1*kCexh + d2*i3 + d3*i4 + i5] *
            ((data1[d1*kHz + d2*i3 + d3*i4 + i5]
              - data1[d1*kHz + d2*i3 + d3*(i4-1) + i5])
             - (data1[d1*kHy + d2*i3 + d3*i4 + i5]
                - data1[d1*kHy + d2*i3 + d3*i4 + i5 - 1]));

        data1[d1*kEx + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kCexe + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kEx + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kCexh + d2*i3 + d3*i4 + i5+1] *
            ((data1[d1*kHz + d2*i3 + d3*i4 + i5+1]
              - data1[d1*kHz + d2*i3 + d3*(i4-1) + i5+1])
	   - (data1[d1*kHy + d2*i3 + d3*i4 + i5+1]
              - data1[d1*kHy + d2*i3 + d3*i4 + i5+1 - 1]));
      }

      for (int i5 = z1; i5 <= z2; i5+=2) {
        data1[d1*kEy + d2*i3 + d3*i4 + i5]
            = data0[d1*kCeye + d2*i3 + d3*i4 + i5] *
            data0[d1*kEy + d2*i3 + d3*i4 + i5]
            + data0[d1*kCeyh + d2*i3 + d3*i4 + i5] *
            ((data1[d1*kHx + d2*i3 + d3*i4 + i5]
              - data1[d1*kHx + d2*i3 + d3*i4 + i5 - 1])
             - (data1[d1*kHz + d2*i3 + d3*i4 + i5]
                - data1[d1*kHz + d2*(i3-1) + d3*i4 + i5]));

        data1[d1*kEy + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kCeye + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kEy + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kCeyh + d2*i3 + d3*i4 + i5+1] *
            ((data1[d1*kHx + d2*i3 + d3*i4 + i5+1]
              - data1[d1*kHx + d2*i3 + d3*i4 + i5+1 - 1])
	   - (data1[d1*kHz + d2*i3 + d3*i4 + i5+1]
              - data1[d1*kHz + d2*(i3-1) + d3*i4 + i5+1]));
      }

      for (int i5 = z1; i5 <= z2; i5+=2) {
	data1[d1*kEz + d2*i3 + d3*i4 + i5]
            = data0[d1*kCeze + d2*i3 + d3*i4 + i5] *
            data0[d1*kEz + d2*i3 + d3*i4 + i5]
            + data0[d1*kCezh + d2*i3 + d3*i4 + i5] *
            ((data1[d1*kHy + d2*i3 + d3*i4 + i5]
              - data1[d1*kHy + d2*(i3-1) + d3*i4 + i5])
             - (data1[d1*kHx + d2*i3 + d3*i4 + i5]
                - data1[d1*kHx + d2*i3 + d3*(i4-1) + i5]))
            + data0[d1*kSrcEz + d2*i3 + d3*i4 + i5];
        ++num_of_updates;

        data1[d1*kEz + d2*i3 + d3*i4 + i5+1]
            = data0[d1*kCeze + d2*i3 + d3*i4 + i5+1] *
            data0[d1*kEz + d2*i3 + d3*i4 + i5+1]
            + data0[d1*kCezh + d2*i3 + d3*i4 + i5+1] *
	  ((data1[d1*kHy + d2*i3 + d3*i4 + i5+1]
            - data1[d1*kHy + d2*(i3-1) + d3*i4 + i5+1])
	   - (data1[d1*kHx + d2*i3 + d3*i4 + i5+1]
              - data1[d1*kHx + d2*i3 + d3*(i4-1) + i5+1]))
            + data0[d1*kSrcEz + d2*i3 + d3*i4 + i5+1];
       ++num_of_updates;
      }            
    }
  }
}
void InitializeFields() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z + 1; ++i5) {
        data0[d1*kHx + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kHy + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kHz + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kEx + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kEy + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kEz + d2*i3 + d3*i4 + i5] = initial;
        data0[d1*kSrcEz + d2*i3 + d3*i4 + i5] = initial;

        data1[d1*kHx + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kHy + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kHz + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kEx + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kEy + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kEz + d2*i3 + d3*i4 + i5] = initial;
        data1[d1*kSrcEz + d2*i3 + d3*i4 + i5] = initial;
      }
    }
  }
}

void SetCoefficients() {
  for (int i3 = 0; i3 <= length_x + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z + 1; ++i5) {
        data0[d1*kChxh + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kChxe + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kChyh + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kChye + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kChzh + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kChze + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kCexe + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kCexh + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kCeye + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kCeyh + d2*i3 + d3*i4 + i5] = coefficient;
        data0[d1*kCeze + d2*i3 + d3*i4 + i5] = 1.0;
        data0[d1*kCezh + d2*i3 + d3*i4 + i5] = coefficient;

        data1[d1*kChxh + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kChxe + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kChyh + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kChye + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kChzh + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kChze + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kCexe + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kCexh + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kCeye + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kCeyh + d2*i3 + d3*i4 + i5] = coefficient;
        data1[d1*kCeze + d2*i3 + d3*i4 + i5] = 1.0;
        data1[d1*kCezh + d2*i3 + d3*i4 + i5] = coefficient;
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
	sum += ABS(data0[d1*kEz + d2*i3 + d3*i4 + i5]);
      }
    }
  }
  return sum;
}

double UsingStdC() {
  SetCoefficients();
  InitializeFields();
  double time = 0, start, end;
  start = MPI_Wtime();
  clock_t start_time = clock();
  double src;
  for (long int t = 0; t < steps; ++t) {
    src = exp(-1.0 * (double) (t) * (double) (t) / ((double) steps * (double) steps));
    for (int i3 = 0; i3 <= length_x + 1; ++i3) {
      for (int i4 = 0; i4 <= length_y + 1; ++i4) {
        for (int i5 = 0; i5 <= length_z + 1; ++i5) {
          data0[d1*kSrcEz + d2*i3 + d3*i4 + i5] = src;
        }
      }
    }
    UpdateFields3DStdC(1, length_x, 1, length_y, 1, length_z);
    CycleArrays();
  }
  clock_t finish_time = clock();
  end = MPI_Wtime();
  time += end - start;

  double total_time = ((double)(finish_time - start_time))/CLOCKS_PER_SEC;
  if (time_flag == 1) {
    printf(">> 3D StdC %dx%dx%dx%ld(steps) took %f \n",
           length_x, length_y, length_z, steps, total_time);
  }
  if (sum_flag == 1) {
    printf(">> 3D StdC %dx%dx%dx%ld(steps) control sum for Ez is %23.16e\n",
           length_x, length_y, length_z, steps, EzControlSum());
  }
  if (size_check_flag == 1)
   printf("================= Updates: %li of %li (ratio %10.8f) ==============\n",
          num_of_updates, (SIZE*SIZE*SIZE*STEPS),
          (double)num_of_updates/(double)(SIZE*SIZE*SIZE*STEPS));
  return total_time;
}

void RunCalculation() {
  initFDTD();
  UsingStdC();
  freeFDTD();
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int total_processes_number, process_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  RunCalculation();
  MPI_Finalize();
  return 0;
}

