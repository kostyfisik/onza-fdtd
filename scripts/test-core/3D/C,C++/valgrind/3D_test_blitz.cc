/// @file   3D_test_blitz.cc
/// @author Markovich Dmitry <dmmrkovich at gmail (.) com>
/// @copyright 2012 Markovich Dmitry
#include <blitz/array.h>
#include <blitz/array/stencilops.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#define SM (CLS / sizeof(double))

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

const int depth = 2;
int components, repeats, min_size, max_size, max_steps;
int sum_flag, print_flag, time_flag, block_flag;
blitz::Array <blitz::Array <double, 4>, 1> data(depth);
blitz::Range i = blitz::Range::all(), j = blitz::Range::all(),
  k = blitz::Range::all();
blitz::Array<double, 1> time_stencil;
blitz::Array<double, 2> time_stencil_blocks;
blitz::Array<double, 1> time_range;
blitz::Array<double, 2> time_range_blocks;
int length_x, length_y, length_z, steps;
double initial = 0.0;
enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze,
                     kSrcEz };

int CheckCLS() {
  /*
  printf("\n>> Default CLS(cache line size) is %d (bytes)\n", CLS);
  printf("\n>> Default SM(cache line size) is %d (elements)\n", SM);
  */
  int type;
  {
    double dummy = 0.0;
    type = sizeof(dummy);
  }
  int cls = CLS, sm = SM;
  if (cls == 0) {
    cls = 64;
    sm = cls / type;
    /*
    printf("\n>> Default CLS is 0. Will automatically set values.\n");
    printf("\n>> CLS is set to %d (bytes)\n", cls);
    printf("\n>> SM is set to %d (elements)\n", sm);
    */
  }
  return sm;
}

void UpdateFields3DStencil(blitz::Range i, blitz::Range j, blitz::Range k) {
 
  {
    const blitz::Range j1(j(0), 1 + j(j.length() - 1)),
      k1(k(0), 1 + k(k.length() - 1));
    blitz::Array<double, 3> tmpkHx_next = (data(1)(kHx, i, j1, k1));
    blitz::Array<double, 3> tmpkChxh = (data(0)(kChxh, i, j1, k1));
    blitz::Array<double, 3> tmpkHx = (data(0)(kHx, i, j1, k1));
    blitz::Array<double, 3> tmpkChxe = (data(0)(kChxe, i, j1, k1));
    blitz::Array<double, 3> tmpkEy = (data(0)(kEy, i, j1, k1));
    blitz::Array<double, 3> tmpkEz = (data(0)(kEz, i, j1, k1));
    applyStencil(update_field_Hx(), tmpkHx_next, tmpkChxh, tmpkHx,
                 tmpkChxe, tmpkEy, tmpkEz);
    data(1)(kHx, i, j1, k1) = tmpkHx_next;
  }
  {
    const blitz::Range i1(i(0), 1 + i(i.length() - 1)),
      k1(k(0), 1 + k(k.length() - 1));
    blitz::Array<double, 3> tmpkHy_next = (data(1)(kHy, i1, j, k1));
    blitz::Array<double, 3> tmpkChyh = (data(0)(kChyh, i1, j, k1));
    blitz::Array<double, 3> tmpkHy = (data(0)(kHy, i1, j, k1));
    blitz::Array<double, 3> tmpkChye = (data(0)(kChye, i1, j, k1));
    blitz::Array<double, 3> tmpkEz = (data(0)(kEz, i1, j, k1));
    blitz::Array<double, 3> tmpkEx = (data(0)(kEx, i1, j, k1));
    applyStencil(update_field_Hy(), tmpkHy_next, tmpkChyh, tmpkHy,
                 tmpkChye, tmpkEz, tmpkEx);
    data(1)(kHy, i1, j, k1) = tmpkHy_next;
  }
  {
  const blitz::Range i1(i(0), 1 + i(i.length() - 1)),
      j1(j(0), 1 + j(j.length() - 1));
    blitz::Array<double, 3> tmpkHz_next = (data(1)(kHz, i1, j1, k));
    blitz::Array<double, 3> tmpkChzh = (data(0)(kChzh, i1, j1, k));
    blitz::Array<double, 3> tmpkHz = (data(0)(kHz, i1, j1, k));
    blitz::Array<double, 3> tmpkChze = (data(0)(kChze, i1, j1, k));
    blitz::Array<double, 3> tmpkEx = (data(0)(kEx, i1, j1, k));
    blitz::Array<double, 3> tmpkEy = (data(0)(kEy, i1, j1, k));
    applyStencil(update_field_Hz(), tmpkHz_next, tmpkChzh, tmpkHz,
                 tmpkChze, tmpkEx, tmpkEy);
    data(1)(kHz, i1, j1, k) = tmpkHz_next;
  }
  {
    const blitz::Range j1(-1 + j(0), j(j.length() - 1)),
      k1(-1 + k(0), k(k.length() - 1));
    blitz::Array<double, 3> tmpkEx_next = (data(1)(kEx, i, j1, k1));
    blitz::Array<double, 3> tmpkCexe = (data(0)(kCexe, i, j1, k1));
    blitz::Array<double, 3> tmpkEx = (data(0)(kEx, i, j1, k1));
    blitz::Array<double, 3> tmpkCexh = (data(0)(kCexh, i, j1, k1));
    blitz::Array<double, 3> tmpkHz = (data(0)(kHz, i, j1, k1));
    blitz::Array<double, 3> tmpkHy = (data(0)(kHy, i, j1, k1));
    applyStencil(update_field_Ex(), tmpkEx_next, tmpkCexe, tmpkEx,
                 tmpkCexh, tmpkHz, tmpkHy);
    data(1)(kEx, i, j1, k1) = tmpkEx_next;
  }
  {
    const blitz::Range i1(-1 + i(0), i(i.length() - 1)),
      k1(-1 + k(0), k(k.length() - 1));
    blitz::Array<double, 3> tmpkEy_next = (data(1)(kEy, i1, j, k1));
    blitz::Array<double, 3> tmpkCeye = (data(0)(kCeye, i1, j, k1));
    blitz::Array<double, 3> tmpkEy = (data(0)(kEy, i1, j, k1));
    blitz::Array<double, 3> tmpkCeyh = (data(0)(kCeyh, i1, j, k1)); 
    blitz::Array<double, 3> tmpkHx = (data(0)(kHx, i1, j, k1));
    blitz::Array<double, 3> tmpkHz = (data(0)(kHz, i1, j, k1));
    applyStencil(update_field_Ey(), tmpkEy_next, tmpkCeye, tmpkEy,
                 tmpkCeyh, tmpkHx, tmpkHz);
    data(1)(kEy, i1, j, k1) = tmpkEy_next;
  }
  {
    const blitz::Range i1(-1 + i(0), i(i.length() - 1)),
      j1(-1 + j(0), j(j.length() - 1));  
    blitz::Array<double, 3> tmpkEz_next = (data(1)(kEz, i1, j1, k));
    blitz::Array<double, 3> tmpkCeze = (data(0)(kCeze, i1, j1, k));
    blitz::Array<double, 3> tmpkEz = (data(0)(kEz, i1, j1, k));
    blitz::Array<double, 3> tmpkCezh = (data(0)(kCezh, i1, j1, k));
    blitz::Array<double, 3> tmpkSrcEz = (data(0)(kSrcEz, i1, j1, k));
    blitz::Array<double, 3> tmpkHy = (data(0)(kHy, i1, j1, k));
    blitz::Array<double, 3> tmpkHx = (data(0)(kHx, i1, j1, k));
    applyStencil(update_field_Ez(), tmpkEz_next, tmpkCeze, tmpkEz,
                 tmpkCezh, tmpkHy, tmpkHx, tmpkSrcEz);
    data(1)(kEz, i1, j1, k) = tmpkEz_next;
  }
}
double UsingStencils(int length_x, int length_y, int length_z,
                     int steps) {
  double time = 0, start, end;
  const blitz::Range first(1, length_x), second(1, length_y),
    third(1, length_z);
  for (int o = 0; o < repeats; ++o) {
    data(0)(kHx, i, j, k) = initial;
    data(0)(kHy, i, j, k) = initial;
    data(0)(kHz, i, j, k) = initial;
    data(0)(kEx, i, j, k) = initial;
    data(0)(kEy, i, j, k) = initial;
    data(0)(kEz, i, j, k) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src = exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j, k) = src;
      UpdateFields3DStencil(first, second, third);
      cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Stencil 3D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Stencil 3D control sum for Ez is "
              << sum(abs(data(0)(static_cast<int>(kEz), i, j, k))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
}

double RunAllBlocksForModelStencil(int length_x, int length_y, int length_z,
				   int max, int steps, int r) {
  const int lim = ceil(max / r);
  const int sm_x_ = length_x / lim, sm_y_ = length_y / lim,
    sm_z_ = length_z / lim;
  double time = 0, start, end;
  time = 0;
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
    data(0)(kHx, i, j, k) = initial;
    data(0)(kHy, i, j, k) = initial;
    data(0)(kHz, i, j, k) = initial;
    data(0)(kEx, i, j, k) = initial;
    data(0)(kEy, i, j, k) = initial;
    data(0)(kEz, i, j, k) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j, k) = src;
      for (int v1 = 1; v1 <= lim; ++v1) {
	for (int v2 = 1; v2 <= lim; ++v2) {
	  for (int v3 = 1; v3 <= lim; ++v3) {
	    UpdateFields3DStencil(Block_ranges_x(v1 - 1),
				  Block_ranges_y(v2 - 1),
				  Block_ranges_z(v3 - 1));
	  }
	}
      }
      cycleArrays(data(0), data(1));
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
    std::cout << sum(abs(data(0)(static_cast<int>(kEz), i, j, k))) << "\n";
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
}

void UpdateFields3DRange(blitz::Range i, blitz::Range j, blitz::Range k) {
  data(1)(kHx, i, j, k) = data(0)(kChxh, i, j, k) *
          data(0)(kHx, i, j, k) + data(0)(kChxe, i, j, k) *
          ((data(0)(kEy, i, j, k + 1) - data(0)(kEy, i, j, k))
           - (data(0)(kEz, i, j + 1, k) - data(0)(kEz, i, j, k)));
 
  data(1)(kHy, i, j, k) = data(0)(kChyh, i, j, k) *
          data(0)(kHy, i, j, k) + data(0)(kChye, i, j, k) *
          ((data(0)(kEz, i + 1, j, k) - data(0)(kEz, i, j, k))
           - (data(0)(kEx, i, j, k + 1) - data(0)(kEx, i, j, k)));
  
  data(1)(kHz, i, j, k) = data(0)(kChzh, i, j, k) *
          data(0)(kHz, i, j, k) + data(0)(kChze, i, j, k) *
          ((data(0)(kEx, i, j + 1, k) - data(0)(kEx, i, j, k))
           - (data(0)(kEy, i + 1, j, k) - data(0)(kEy, i, j, k)));
  
  data(1)(kEx, i, j, k) = data(0)(kCexe, i, j, k) *
          data(0)(kEx, i, j, k) + data(0)(kCexh, i, j, k) *
          ((data(0)(kHz, i, j, k) - data(0)(kHz, i, j - 1, k))
           - (data(0)(kHy, i, j, k) - data(0)(kHy, i, j, k - 1)));
  
  data(1)(kEy, i, j, k) = data(0)(kCeye, i, j, k) *
          data(0)(kEy, i, j, k) + data(0)(kCeyh, i, j, k) *
          ((data(0)(kHx, i, j, k) - data(0)(kHx, i, j, k - 1))
           - (data(0)(kHz, i, j, k) - data(0)(kHy, i - 1, j, k)));

  data(1)(kEz, i, j, k) = data(0)(kCeze, i, j, k) *
          data(0)(kEz, i, j, k) + data(0)(kCezh, i, j, k) *
          ((data(0)(kHy, i, j,  k) - data(0)(kHy, i - 1, j, k))
           - (data(0)(kHx, i, j, k) - data(0)(kHx, i, j - 1, k)))
          + data(0)(kSrcEz, i, j, k);
}

double UsingRange(int length_x, int length_y, int length_z, int steps) {
  const blitz::Range inner_i(1, length_x);
  const blitz::Range inner_j(1, length_y);
  const blitz::Range inner_k(1, length_z);
  double time = 0, start, end;
  for (int l = 0; l < repeats; ++l) {
    data(0)(kHx, i, j, k) = initial;
    data(0)(kHy, i, j, k) = initial;
    data(0)(kHz, i, j, k) = initial;
    data(0)(kEx, i, j, k) = initial;
    data(0)(kEy, i, j, k) = initial;
    data(0)(kEz, i, j, k) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =  exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j, k) = src;
      UpdateFields3DRange(inner_i, inner_j, inner_k);
            cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Range 3D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Range 3D control sum for Ez is "
              << sum(abs(data(0)(static_cast<int>(kEz), i, j, k))) << "\n";
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
}

double RunAllBlocksForModelRange(int length_x, int length_y, int length_z,
                                 int max, int steps, int r) {
  int lim = ceil(max / r);
  const int sm_x_ = length_x / lim, sm_y_ = length_y / lim,
    sm_z_ = length_z / lim;
  double time = 0, start, end;
  time = 0;
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
    data(0)(kHx, i, j, k) = initial;
    data(0)(kHy, i, j, k) = initial;
    data(0)(kHz, i, j, k) = initial;
    data(0)(kEx, i, j, k) = initial;
    data(0)(kEy, i, j, k) = initial;
    data(0)(kEz, i, j, k) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j, k) = src;
      for (int v1 = 1; v1 <= lim; ++v1) {
       for (int v2 = 1; v2 <= lim; ++v2) {
         for (int v3 = 1; v3 <= lim; ++v3) {
           UpdateFields3DRange(Block_ranges_x(v1 - 1),
			       Block_ranges_y(v2 - 1),
                               Block_ranges_z(v3 - 1));
         }
       }
      }
      cycleArrays(data(0), data(1));
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
    std::cout << sum(abs(data(0)(static_cast<int>(kEz), i, j, k))) << "\n";
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kEz), i, j, k);
  }
  return (time / repeats);
}

void PreparePlotDataFile(int length_x, int length_y, int length_z, int max,
                         int steps) {
  char basename[80] = "blitz_size_", filename[100];
  if (max < 100) {
    snprintf(filename, sizeof(filename), "%s00%dx00%dx00%d.plot",
             basename, length_x, length_y, length_z);
  } else {
    if (max < 1000) {
      snprintf(filename, sizeof(filename), "%s0%dx0%dx0%d.plot",
               basename, length_x, length_y, length_z);
    } else {
      snprintf(filename, sizeof(filename), "%s%dx%dx%d.plot",
               basename, length_x, length_y, length_z);
    }
  }
  FILE *f;
  f = fopen(filename, "a");
  fprintf(f, "#Model dimensions: %dx%dx%dx%d(steps)\n",
          length_x, length_y, length_z, steps);
  if (block_flag == 0) {
    fprintf(f, "#Blocks' length\tStencil");
    fprintf(f, "\tRange\t");
    fprintf(f, "\t\tt(stncl)\t\tt(rng)\n");
  }
  else {
    fprintf(f, "#Blocks' length\tStencil\tStencil(blocks)");
    fprintf(f, "\tRange\tRange(blocks)\tBest");
    fprintf(f, "\t\tt(stncl)\tt(stncl_blcks)\tt(rng)\t\tt(rng_blcks)\n");
  }
  fclose(f);
}

void WritePlotData(int length_x, int length_y, int length_z, int max,
                   int blocks, double t1, double t2, double t3, double t4,
                   int steps) {
  char basename[80] = "blitz_size_", filename[100];
  if (max < 100) {
    snprintf(filename, sizeof(filename), "%s00%dx00%dx00%d.plot",
             basename, length_x, length_y, length_z);
  } else {
    if (max < 1000) {
      snprintf(filename, sizeof(filename), "%s0%dx0%dx0%d.plot",
               basename, length_x, length_y, length_z);
    } else {
      snprintf(filename, sizeof(filename), "%s%dx%dx%d.plot",
               basename, length_x, length_y, length_z);
    }
  }
  FILE *f;
  f = fopen(filename, "a");
  if (block_flag == 0) {
    int total = length_x * length_y * length_z * steps;
    fprintf(f, "\t%.2f\t%.2f\t", t1, t3);
    fprintf(f, "\t%.2e\t%.2e\n", (t1 / total), (t3 / total));
  }
  else {
  fprintf(f, "%d\t\t%.2f\t%.2f\t\t%.2f\t%.2f\t", blocks, t1, t3, t4);
  double min = t1;
  min = (t2 < min) ? t2:min;
  min = (t3 < min) ? t3:min;
  min = (t4 < min) ? t4:min;
  if (min == t1) {
    fprintf(f, "\tStencil(+%.1f%)\t", (100 * (t3-t1) / t3));
  }
  if (min == t2) {
    fprintf(f, "\tStencil(blocks)(+%.1f%)\t", (100 * (t3-t2) / t3));
  }
  if (min == t3) {
    fprintf(f, "\tRange\t");
  }
  if (min == t4) {
    fprintf(f, "\tRange(blocks)(+%.1f%)\t", (100 * (t3-t4) / t3));
  }
  int total = length_x * length_y * length_z * steps;
  fprintf(f, "\t%.2e\t%.2e\t%.2e\t%.2e\n", (t1 / total), (t2 / total),
          (t3 / total), (t4 / total));
  }
  fclose(f);
}
int GetPowerOfTwo(int number) {
  int result = number;
  int power = 0;
  while (result != 1) {
    result = result / 2;
    power += 1;
  }
  return power;
}
void RunAllModelsCalculation() {
  int sm_min = CheckCLS(), sm_max = max_size;
  int C = max_steps * max_size * max_size * max_size;
  int counter1 = 0, counter2;
  time_stencil.resize(1 + GetPowerOfTwo(max_size / min_size));
  time_stencil_blocks.resize((1 + GetPowerOfTwo(max_size / min_size)),
                           (GetPowerOfTwo(sm_max / sm_min)));
  time_range.resize(1 + GetPowerOfTwo(max_size / min_size));
  time_range_blocks.resize((1 + GetPowerOfTwo(max_size / min_size)),
                           (GetPowerOfTwo(sm_max / sm_min)));
  blitz::Array <int, 1> blocks(GetPowerOfTwo(sm_max / sm_min));
  for (int size = min_size; size <= max_size; size *=2) {
    length_x = size;
    length_y = size;
    length_z = size;
    steps = C / (length_x * length_y * length_z);
    int max = (length_x <= length_y) ? length_y:length_x;
    max = (length_z <= max) ? max:length_z;
    PreparePlotDataFile(length_x, length_y, length_z, max, steps);
    for (int counter = 0; counter < depth; ++counter) {
      data(counter).resize(components, length_x + 2, length_y + 2,
                           length_z + 2);
      data(counter) = 0;
    }
    double coefficient = 0.1;
    for (int w = 0; w < depth; w++) {
      data(w)(kCexe, i, j, k) = coefficient;
      data(w)(kCexh, i, j, k) = coefficient;
      data(w)(kCeye, i, j, k) = coefficient;
      data(w)(kCeyh, i, j, k) = coefficient;
      data(w)(kCeze, i, j, k) = 1;
      data(w)(kCezh, i, j, k) = coefficient;
      data(w)(kChxe, i, j, k) = coefficient;
      data(w)(kChxh, i, j, k) = coefficient;
      data(w)(kChye, i, j, k) = coefficient;
      data(w)(kChyh, i, j, k) = coefficient;
      data(w)(kChze, i, j, k) = coefficient;
      data(w)(kChzh, i, j, k) = coefficient;
    }
      printf("\n>> __Processing task: 3D %dx%dx%dx%d(steps)\n",
             length_x, length_y, length_z, steps);
    time_stencil(counter1) = UsingStencils(length_x, length_y,
					   length_z, steps);
    time_range(counter1) = UsingRange(length_x, length_y, length_z, steps);
    
    if (sm_min < max) {
      counter2 = 0;
      for (int r = sm_min; r <= max; r *= 2) {
        blocks(counter2) = r;
        if (r < max) {
          if (block_flag == 0) {
            time_stencil_blocks(counter1, counter2) = 0.0;
            time_range_blocks(counter1, counter2) = 0.0;
          }
          else {
            time_stencil_blocks(counter1, counter2) =
              RunAllBlocksForModelStencil(length_x, length_y, length_z, max,
                                          steps, r);
            time_range_blocks(counter1, counter2) =
              RunAllBlocksForModelRange(length_x, length_y, length_z, max,
                                        steps, r);
          }
        }
        counter2 += 1;
      }
      for (int cnt = 0; cnt < (counter2 - 1); ++cnt) {
	WritePlotData(length_x, length_y, length_z, max, blocks(cnt),
		      time_stencil(counter1),
		      time_stencil_blocks(counter1, cnt),
		      time_range(counter1),
		      time_range_blocks(counter1, cnt), steps);
      }
    }  // end of if (sm_min < max)
    counter1 += 1;
      printf("\n>> __Done processing task: 3D %dx%dx%dx%d(steps)\n",
             length_x, length_y, length_z, steps);
  }  // end of for(int size = min_size; size <= max_size; size *= 2)
}
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int total_processes_number, process_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  // Fixed parameters **********************************************************
  components = kSrcEz + 1;
  sum_flag = 0;
  print_flag = 0;
  time_flag = 1;
  block_flag = blk;
  int size_check_flag = 0;
  // ***************************************************************************
  // Parameters to set *********************************************************
  // set number of timesteps and number of time measure procedure repeats
  max_steps = 10;
  repeats = 2;
  // set min ( >= 16 ) and max size of model
  min_size = 32;
  max_size = 128;
  // ***************************************************************************
  if (size_check_flag == 1) {
    if (min_size < 16) {
      printf("\n>> Automatically setting min_size = 16\n");
      min_size = 16;
    }
    if (max_size < min_size) {
      printf("\n>> Automatically setting max_size = 32\n");
      max_size = 32;
    }
  }
  RunAllModelsCalculation();
  MPI_Finalize();
  return 0;
}

