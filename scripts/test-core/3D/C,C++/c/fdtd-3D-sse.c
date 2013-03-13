/// @file   3D_test_std.cc
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
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <xmmintrin.h>

#define SIZE 128
#define STEPS ((int)128*128*128*100/SIZE/SIZE/SIZE)
/* #define SIZE 4 */
/* #define STEPS 4 */
#define full_output 0
// vec - 2 for sse, 4 for avx
#define vec 2
const int size = SIZE;
int steps = STEPS;

const int depth = 2;
const int length_x = SIZE;
const int length_y = SIZE;
const int length_z = SIZE;

const int full_output_flag = full_output;
int sum_flag, time_flag, size_check_flag;
const double initial = 0.000, boundary = 0.0, coefficient = 0.1;
enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze, 
                     kSrcEz };
char* name[]= {"kEx", "kCexe", "kCexh",
    "kEy", "kCeye", "kCeyh",
    "kEz", "kCeze", "kCezh",
    "kHx", "kChxh", "kChxe",
    "kHy", "kChyh", "kChye",
    "kHz", "kChzh", "kChze",
    "kSrcEz" };
const int number_of_components = kSrcEz + 1;
const int s_z = 1;
const int s_y = SIZE + vec*2;
const int s_x = (SIZE + vec*2)*(SIZE + vec*2);
const int s_c = (SIZE + vec*2)*(SIZE + vec*2)*(SIZE + vec*2);
/* double *__restrict data0; */
/* double *__restrict data1; */
void UpdateFields3DStdC();
void InitializeFields();
void SetCoefficients();
void CycleArrays();
double EzControlSum();
double UsingStdC();

//debug
const int k = SIZE+2;
/* double *x; double *y; // vectors of length k */
/***************************************************************************/  
/****************   Memory managment  **************************************/  
/***************************************************************************/  
struct AlignedMem
{
  int aligment_;  // in bytes, usually 16 (sse), 32 (avx) or 64 (xeon phi SIMD)
  // double is 8 bytes, so sse register can fit 2 double numbers.
  void * unaligned_ptr_;
  // TODO (check if restricted is needed)
  double *__restrict d_; // aligned pointer
  //double * d_; // aligned pointer
} data0, data1, x, y;
int AlignedMalloc(struct AlignedMem* data) {
  /* After answer of Jonathan Leffler in http://stackoverflow
     .com/questions/227897
     /solve-the-memory-alignment-in-c-interview-question-that-stumped-me
     With specialization to our FDTD memory model.
  */
  /* int aligment_in_bytes = 16; */
  /* int aligment_in_bytes = 32; */
  int aligment_in_bytes = 8*vec;
  data->aligment_ = aligment_in_bytes;
  if (sizeof(double)!= 8) {
    printf("Unsupported size of double %i (should be 8)\n", sizeof(double));
    exit(1);
  }
  data->unaligned_ptr_ = malloc(sizeof(double)*number_of_components*s_c
                                                 + data->aligment_);
  uintptr_t mask = ~(uintptr_t)(data->aligment_-1);
  /* printf("mask: %p\n", mask); */
  data->d_ = (double *)(((uintptr_t)data->unaligned_ptr_
                         + data->aligment_) & mask);
  return 1;
}
int AlignedFree(struct AlignedMem* data) {
  free(data->unaligned_ptr_);
  return 1;
}
void InitFDTD() {
  AlignedMalloc(&data0); 
  AlignedMalloc(&data1);
  //debug
  AlignedMalloc(&x);
  AlignedMalloc(&y);

  /* printf("data0.u_ addr: %p\n",data0.unaligned_ptr_);  */
  /* printf("data0.d_ addr: %p\n",data0.d_);  */
}
void  FreeFDTD() {
  AlignedFree(&data0);
  AlignedFree(&data1);
  //debug
  AlignedFree(&x);
  AlignedFree(&y);
}
/***************************************************************************/  
/****************  End of memory managment block  **************************/  
/***************************************************************************/  

/***************************************************************************/  
/***************************************************************************/  
/***************************************************************************/  

double ABS(double number) {
  if (number >= 0) {
    return number;
  } else {
    return -1 * number;
  }
}
// Prefetch definition taken from:
// http://software.intel.com/en-us/forums/showthread.php?t=46284
// tune FETCH_DISTANCE according to real world experiments
#define PREFETCH_T0(addr, nrOfBytesAhead) _mm_prefetch(((char *)(addr))+nrOfBytesAhead, _MM_HINT_T0)
#define FETCH_DISTANCE        256*2
/* PREFETCH_T0( &data0.d_[s_c*kChxh + s_x*i_x + s_y*i_y     + 2], FETCH_DISTANCE); */

void UpdateFields3DStdC() {
  __m128d m1, m2, m3,m4;
  __m128d t1kHx, t0kHx, t1kHy, t0kHy, t1kHz, t0kHz,
      t1kEx, t0kEx, t1kEy, t0kEy, t1kEz, t0kEz,
      t0kChxh, t0kChxe, t0kChyh, t0kChye, t0kChzh, t0kChze,
      t0kCexh, t0kCexe, t0kCeyh, t0kCeye, t0kCezh, t0kCeze,
      t0kEx1py, t0kEx1pz, t0kEy1px, t0kEy1pz, t0kEz1px, t0kEz1py,
      t1kHx1my, t1kHx1mz, t1kHy1mx, t1kHy1mz, t1kHz1mx, t1kHz1my,
      t0kSrcEz;
  /* __m128d *pt1kHx, *pt0kHx, *pt1kHy, *pt0kHy, *pt1kHz, *pt0kHz, */
  /*     *pt1kEx, *pt0kEx, *pt1kEy, *pt0kEy, *pt1ykEz, *pt0kEz, */
  /*     *pt0kChxh, *pt0kChxe, *pt0kChyh, *pt0kChye, *pt0kChzh, *pt0kChze, */
  /*     *pt0kCexh, *pt0kCexe, *pt0kCeyh, *pt0kCeye, *pt0kCezh, *pt0kCeze, */
  /*     *pt0kEx1py, *pt0kEx1pz, *pt0kEy1px, *pt0kEy1pz, *pt0kEz1px, *pt0kEz1py, */
  /*     *pt1kHx1my, *pt1kHx1mz, *pt1kHy1mx, *pt1kHy1mz, *pt1kHz1mx, *pt1kHz1my, */
  /*     *pt0kSrcEz; */
  for (int i_x = vec; i_x < length_x+vec; ++i_x) {
    for (int i_y = vec; i_y < length_y+vec; ++i_y) {
      const int cell_z0 = s_x*i_x + s_y*i_y;
      for(int i_z=vec;i_z<length_z+vec;i_z+=vec){        
        t0kEy1pz  = _mm_loadu_pd(&data0.d_[s_c*kEy   + cell_z0 + s_z+i_z]);
        t0kChxh  = _mm_load_pd (&data0.d_[s_c*kChxh + cell_z0      +i_z]);
        t0kHx    = _mm_load_pd (&data0.d_[s_c*kHx   + cell_z0      +i_z]);
        t0kChxe  = _mm_load_pd (&data0.d_[s_c*kChxe + cell_z0      +i_z]);
        t0kEy    = _mm_load_pd (&data0.d_[s_c*kEy   + cell_z0      +i_z]);
        t0kEz1py = _mm_load_pd (&data0.d_[s_c*kEz   + cell_z0 + s_y+i_z]);
        t0kEz    = _mm_load_pd (&data0.d_[s_c*kEz   + cell_z0      +i_z]);
        t1kHx    = _mm_load_pd (&data1.d_[s_c*kHx   + cell_z0      +i_z]);
        t1kHx=t0kChxh*t0kHx+t0kChxe*((t0kEy1pz-t0kEy)-(t0kEz1py-t0kEz));
        _mm_store_pd(&data1.d_[s_c*kHx+cell_z0+i_z],t1kHx);

        /*t1kHy=t0kChyh*t0kHy+t0kChye*((t0kEz1px-t0kEz)-(t0kEx1pz-t0kEx));*/
        t0kChyh=_mm_load_pd(&data0.d_[s_c*kChyh+cell_z0+i_z]);
        t0kHy=_mm_load_pd(&data0.d_[s_c*kHy+cell_z0+i_z]);
        t0kChye=_mm_load_pd(&data0.d_[s_c*kChye+cell_z0+i_z]);
        t0kEz1px=_mm_load_pd(&data0.d_[s_c*kEz+s_x*(i_x+1)+s_y*i_y+i_z]);
        /* t0kEz=_mm_load_pd(&data0.d_[s_c*kEz+cell_z0+i_z]); */
        t0kEx1pz=_mm_loadu_pd(&data0.d_[s_c*kEx+cell_z0+i_z+1]);
        t0kEx=_mm_load_pd(&data0.d_[s_c*kEx+cell_z0+i_z]);
        t1kHy=t0kChyh*t0kHy+t0kChye*((t0kEz1px-t0kEz)-(t0kEx1pz-t0kEx));
        _mm_store_pd(&data1.d_[s_c*kHy+cell_z0+i_z],t1kHy);

        /*t1kHz=t0kChzh*t0kHz+t0kChze*((t0kEx1py-t0kEx)-(t0kEy1px-t0kEy));*/
        t0kChzh=_mm_load_pd(&data0.d_[s_c*kChzh+cell_z0+i_z]);
        t0kHz=_mm_load_pd(&data0.d_[s_c*kHz+cell_z0+i_z]);
        t0kChze=_mm_load_pd(&data0.d_[s_c*kChze+cell_z0+i_z]);
        t0kEx1py=_mm_load_pd(&data0.d_[s_c*kEx+s_x*i_x+s_y*(i_y+1)+i_z]);
        /* t0kEx=_mm_load_pd(&data0.d_[s_c*kEx+cell_z0+i_z]); */
        t0kEy1px=_mm_load_pd(&data0.d_[s_c*kEy+s_x*(i_x+1)+s_y*i_y+i_z]);
        /* t0kEy=_mm_load_pd(&data0.d_[s_c*kEy+cell_z0+i_z]); */
        t1kHz=t0kChzh*t0kHz+t0kChze *((t0kEx1py-t0kEx)-(t0kEy1px-t0kEy));
        _mm_store_pd(&data1.d_[s_c*kHz+cell_z0+i_z],t1kHz);
      /* } */
      /* for(int i_z=vec;i_z<length_z+vec;i_z+=vec){ */
        /*t1kEx=t0kCexe*t0kEx+t0kCexh*((t0kHz-t0kHz1my)-(t0kHy-t0kHy1mz));*/
        t0kCexe=_mm_load_pd(&data0.d_[s_c*kCexe+cell_z0+i_z]);
        t0kEx=_mm_load_pd(&data0.d_[s_c*kEx+cell_z0+i_z]);
        t0kCexh=_mm_load_pd(&data0.d_[s_c*kCexh+cell_z0+i_z]);
        /* t1kHz=_mm_load_pd(&data1.d_[s_c*kHz+cell_z0+i_z]); */
        t1kHz1my=_mm_load_pd(&data1.d_[s_c*kHz+s_x*i_x+s_y*(i_y-1)+i_z]);
        /* t1kHy=_mm_load_pd(&data1.d_[s_c*kHy+cell_z0+i_z]); */
        t1kHy1mz=_mm_loadu_pd(&data1.d_[s_c*kHy+cell_z0+i_z-1]);
        t1kEx=t0kCexe*t0kEx+t0kCexh*((t1kHz-t1kHz1my)-(t1kHy-t1kHy1mz));
        _mm_store_pd(&data1.d_[s_c*kEx+cell_z0+i_z],t1kEx);

        /*t1kEy=t0kCeye*t0kEy+t0kCeyh*((t0kHx-t0kHx1mz)-(t0kHz-t0kHz1mx));*/
        t0kCeye  =_mm_load_pd(&data0.d_[s_c*kCeye+cell_z0+i_z]);
        t0kEy    =_mm_load_pd(&data0.d_[s_c*kEy  +cell_z0+i_z]);
        t0kCeyh  =_mm_load_pd(&data0.d_[s_c*kCeyh+cell_z0+i_z]);
        /* t1kHx    =_mm_load_pd(&data1.d_[s_c*kHx  +cell_z0+i_z]); */
        t1kHx1mz =_mm_loadu_pd(&data1.d_[s_c*kHx +cell_z0+i_z-1]);
        /* t1kHz    =_mm_load_pd(&data1.d_[s_c*kHz  +cell_z0+i_z]); */
        t1kHz1mx =_mm_load_pd(&data1.d_[s_c*kHz  +cell_z0 - s_x +i_z]);
        t1kEy = t0kCeye*t0kEy+t0kCeyh*((t1kHx-t1kHx1mz)
                                       -(t1kHz
                                         -t1kHz1mx
                                         )
                                       );
        _mm_store_pd(&data1.d_[s_c*kEy+cell_z0+i_z],t1kEy);

      /*t1kEz=t0kCeze*t0kEz+t0kCezh*((t0kHy-t0kHy1mx)-(t0kHx-t0kHx1my))+t0kSrzEz;*/
        t1kEz=_mm_load_pd(&data1.d_[s_c*kEz+cell_z0+i_z]);
        t0kCeze=_mm_load_pd(&data0.d_[s_c*kCeze+cell_z0+i_z]);
        t0kEz=_mm_load_pd(&data0.d_[s_c*kEz+cell_z0+i_z]);
        t0kCezh=_mm_load_pd(&data0.d_[s_c*kCezh+cell_z0+i_z]);
        /* t1kHy=_mm_load_pd(&data1.d_[s_c*kHy+cell_z0+i_z]); */
        t1kHy1mx=_mm_load_pd(&data1.d_[s_c*kHy+s_x*(i_x-1)+s_y*i_y+i_z]);
        /* t1kHx=_mm_load_pd(&data1.d_[s_c*kHx+cell_z0+i_z]); */
        t1kHx1my=_mm_load_pd(&data1.d_[s_c*kHx+s_x*i_x+s_y*(i_y-1)+i_z]);
        t0kSrcEz=_mm_load_pd(&data0.d_[s_c*kSrcEz+cell_z0+i_z]);
        t1kEz=(t0kCeze*t0kEz+t0kCezh*((t1kHy-t1kHy1mx)-(t1kHx-t1kHx1my)))
              //;
              +t0kSrcEz;
        _mm_store_pd(&data1.d_[s_c*kEz+cell_z0+i_z],t1kEz);
      }
    }
  }
}
void InitializeFields(){
  for(int i_x=0;i_x<length_x+vec*2;++i_x){
    for(int i_y=0;i_y<length_y+vec*2;++i_y){
      for(int i_z=0;i_z<length_z+vec*2;++i_z){
        const int cell = s_x*i_x + s_y*i_y + i_z;
        double set_value = initial;
        if (i_x < vec || i_y < vec || i_z < vec
            || i_x >= length_x + vec
            || i_y >= length_y + vec
            || i_z >= length_z + vec) set_value = boundary;
        data0.d_[s_c*kHx+cell]=set_value;
        data0.d_[s_c*kHy+cell]=set_value;
        data0.d_[s_c*kHz+cell]=set_value;
        data0.d_[s_c*kEx+cell]=set_value;
        data0.d_[s_c*kEy+cell]=set_value;
        data0.d_[s_c*kEz+cell]=set_value;
        //data0.d_[s_c*kSrcEz+cell]=set_value;

        data1.d_[s_c*kHx+cell]=set_value;
        data1.d_[s_c*kHy+cell]=set_value;
        data1.d_[s_c*kHz+cell]=set_value;
        data1.d_[s_c*kEx+cell]=set_value;
        data1.d_[s_c*kEy+cell]=set_value;
        data1.d_[s_c*kEz+cell]=set_value;
        //data1.d_[s_c*kSrcEz+cell]=set_value;
      }
    }
  }
}

void SetCoefficients(){
  for(int i_x=0;i_x<length_x+vec*2;++i_x){
    for(int i_y=0;i_y<length_y+vec*2;++i_y){
      for(int i_z=0;i_z<length_z+vec*2;++i_z){
        const int cell = s_x*i_x + s_y*i_y + i_z;
        data0.d_[s_c*kChxh+cell]=coefficient;
        data0.d_[s_c*kChxe+cell]=coefficient;
        data0.d_[s_c*kChyh+cell]=coefficient;
        data0.d_[s_c*kChye+cell]=coefficient;
        data0.d_[s_c*kChzh+cell]=coefficient;
        data0.d_[s_c*kChze+cell]=coefficient;
        data0.d_[s_c*kCexe+cell]=coefficient;
        data0.d_[s_c*kCexh+cell]=coefficient;
        data0.d_[s_c*kCeye+cell]=coefficient;
        data0.d_[s_c*kCeyh+cell]=coefficient;
        data0.d_[s_c*kCeze+cell]=1;
        data0.d_[s_c*kCezh+cell]=coefficient;

        data1.d_[s_c*kChxh+cell]=coefficient;
        data1.d_[s_c*kChxe+cell]=coefficient;
        data1.d_[s_c*kChyh+cell]=coefficient;
        data1.d_[s_c*kChye+cell]=coefficient;
        data1.d_[s_c*kChzh+cell]=coefficient;
        data1.d_[s_c*kChze+cell]=coefficient;
        data1.d_[s_c*kCexe+cell]=coefficient;
        data1.d_[s_c*kCexh+cell]=coefficient;
        data1.d_[s_c*kCeye+cell]=coefficient;
        data1.d_[s_c*kCeyh+cell]=coefficient;
        data1.d_[s_c*kCeze+cell]=1;
        data1.d_[s_c*kCezh+cell]=coefficient;
      }
    }
  }
}

void CycleArrays(){
  double *tmp;
  tmp=data0.d_;
  data0.d_=data1.d_;
  data1.d_=tmp;
}

double EzControlSum(){
  double sum=0;
  for(int i_x=0;i_x<length_x+vec;++i_x){
    for(int i_y=0;i_y<length_y+vec;++i_y){
      for(int i_z=0;i_z<length_z+vec;++i_z){
	sum+=ABS(data0.d_[s_c*kEz+s_x*i_x+s_y*i_y+i_z]);
      }
    }
  }
  return sum;
}

double UsingStdC(){
  SetCoefficients();
  InitializeFields();
  clock_t start_time=clock();
  /* double src; */
  for(long int t=0;t<steps;++t){
    /* const double src = exp(-(double) (t) * (double) (t) /  */
    /*                        ((double) (steps) * (double) (steps))); */

    const double src =
        exp(-1.0*(double)(t)*(double)(t)/((double)(steps)*(double)(steps))); 
    /* printf("Source = %23.16e \n", src); */

    for(int i_x=vec;i_x<length_x+vec;++i_x){
      for(int i_y=vec;i_y<length_y+vec;++i_y){
        for(int i_z=vec;i_z<length_z+vec;i_z+=1){
          const int cell = s_x*i_x + s_y*i_y + i_z;
          data0.d_[s_c*kSrcEz+cell]=src;
        }
      }
    }
    UpdateFields3DStdC();
    CycleArrays();
  }
  clock_t finish_time=clock();
  double total_time=((double)(finish_time-start_time))/CLOCKS_PER_SEC;
  if(time_flag==1)printf("\n>> 3D sse took%f\n",total_time);
  if (full_output_flag == 1) {
#define format "%23.16e   "
    //#define format "%22.17f "
    //int toPrint = kChze;
    //#define format "%22.13f "
    //int toPrint = kEz;
    //int toPrint = kSrcEz;
    //int toPrint = kEx;
    //int toPrint = kEy;
    //int toPrint = kHz;
    //int toPrint = kHx;
    //int toPrint = kHy;
    int border = vec;
    //int border = 0;
    for (int toPrint = 0; toPrint < number_of_components; ++toPrint) {
      printf("\nSingle sse component %s\n",name[toPrint]);    
      for(int i_x=border;i_x<SIZE+vec*2-border;++i_x){
        for(int i_y=border;i_y<SIZE+vec*2-border;++i_y){
          for(int i_z=border;i_z<SIZE+vec*2-border;i_z+=1){
            const int cell = s_x*i_x + s_y*i_y + i_z;
            printf(format, data0.d_[s_c*toPrint+cell]);
          }
          printf("\n");
        }
        printf("\n\n");
      }
    }
    printf("\n");        
  }
  if(sum_flag==1){
    double control_sum=EzControlSum();
    printf(">> 3D sse %dx%dx%dx%ld(steps) control sum for Ez is %23.16e\n",
           length_x, length_y, length_z, steps, EzControlSum());
  }

  /*printf(">> __Done processing task 3D %dx%dx%dx%d(steps)\n",length_x,length_y,*/
  /*length_z,steps);*/
  return total_time;
}

void RunCalculation(){
  InitFDTD();
  UsingStdC();
  FreeFDTD();
}

int main(int argc,char*argv[]){
  //Flags*********************************************************
  sum_flag=1;
  time_flag=1;
  size_check_flag=1;
  //***************************************************************
  //printf("Mod2\n");
  //debug
  RunCalculation();
  return 0;
}

// StdC 128x128x128x100(steps) control sum
// for Ez is  1.5611249740143090e+08
//            1.5611249321712190e+08
//>> 3D StdC 128x128x128x100(steps) control sum for Ez is  1.5611249740143090e+08
