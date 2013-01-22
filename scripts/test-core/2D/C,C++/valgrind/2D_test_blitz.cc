/// @file   2D_test_blitz.cc
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

const int depth = 2;
int repeats, min_size, max_size, max_steps;
int sum_flag, print_flag, time_flag, block_flag;
blitz::Array <blitz::Array <double, 3>, 1> data(depth);
blitz::Range i = blitz::Range::all(), j = blitz::Range::all();
blitz::Array<double, 1> time_stencil;
blitz::Array<double, 2> time_stencil_blocks;
blitz::Array<double, 1> time_range;
blitz::Array<double, 2> time_range_blocks;
int length_x, length_y, steps;
double initial = 0.0;
enum DataComponents {kEz = 0, kCeze, kCezh,
                     kHy, kChyh, kChye,
                     kHx, kChxh, kChxe, kSrcEz };
int components = kSrcEz + 1;
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
void UpdateFields2DStencil(blitz::Range i, blitz::Range j) {
  {
    const blitz::Range j1(j(0), 1 + j(j.length() - 1)); 
    blitz::Array<double, 2> tmpkHx_next = (data(1)(kHx, i, j1));
    blitz::Array<double, 2> tmpkHx = (data(0)(kHx, i, j1));
    blitz::Array<double, 2> tmpkChxh = (data(0)(kChxh, i, j1));
    blitz::Array<double, 2> tmpkChxe = (data(0)(kChxe, i, j1));
    blitz::Array<double, 2> tmpkEz = (data(0)(kEz, i, j1));
    applyStencil(update_field_Hx(), tmpkHx_next, tmpkHx, tmpkChxh,
                 tmpkChxe, tmpkEz);
    data(1)(kHx, i, j1) = tmpkHx_next;
  }
  {
    const blitz::Range i1(i(0), 1 + i(i.length() - 1));
    blitz::Array<double, 2> tmpkHy_next = (data(1)(kHy, i1, j));
    blitz::Array<double, 2> tmpkHy = (data(0)(kHy, i1, j));
    blitz::Array<double, 2> tmpkChyh = (data(0)(kChyh, i1, j));
    blitz::Array<double, 2> tmpkChye = (data(0)(kChye, i1, j));
    blitz::Array<double, 2> tmpkEz = (data(0)(kEz, i1, j));
    applyStencil(update_field_Hy(), tmpkHy_next, tmpkHy, tmpkChyh,
                 tmpkChye, tmpkEz);
    data(1)(kHy, i1, j) = tmpkHy_next;
  }
  {
    blitz::Array<double, 2> tmpkEz_next = (data(1)(kEz, i, j));
    blitz::Array<double, 2> tmpkEz = (data(0)(kEz, i, j));
    blitz::Array<double, 2> tmpkCeze = (data(0)(kCeze, i, j));
    blitz::Array<double, 2> tmpkCezh = (data(0)(kCezh, i, j));
    blitz::Array<double, 2> tmpkSrcEz = (data(0)(kSrcEz, i, j));
    const blitz::Range j1(-1 + j(0), j(j.length() - 1));
    blitz::Array<double, 2> tmpkHx = (data(0)(kHx, i, j1));
    const blitz::Range i1(-1 + i(0), i(i.length() - 1));
    blitz::Array<double, 2> tmpkHy = (data(0)(kHy, i1, j));
    applyStencil(update_field_Ez(), tmpkEz_next, tmpkEz, tmpkCeze,
                 tmpkCezh, tmpkHy, tmpkHx, tmpkSrcEz);
    data(1)(kEz, i, j) = tmpkEz_next;
  }
}
double UsingStencils(int length_x, int length_y, int steps) {
  double time = 0, start, end;
  const blitz::Range first(1, length_x), second(1, length_y);
  for (int o = 0; o < repeats; ++o) {
    data(0)(kEz, i, j) = initial;
    data(0)(kHx, i, j) = initial;
    data(0)(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j) = src;
      UpdateFields2DStencil(first, second);
      cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Stencil 2D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Stencil 2D control sum for Ez is "
              << sum(abs(data(0)(static_cast<int>(kEz), i, j))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kHx), i, j);
    std::cout << data(0)(static_cast<int>(kHy), i, j);
    std::cout << data(0)(static_cast<int>(kEz), i, j);
  }
  return (time / repeats);
}
double RunAllBlocksForModelStencil(int length_x, int length_y, int steps,
                            int r) {
  int blocks_number_y = ceil(length_y / r);
  int blocks_number_x = ceil(length_x / r);
  int lim, sm_x, sm_y;
  if (blocks_number_x > blocks_number_y) {
    lim = blocks_number_x;
    sm_x = r;
    sm_y = length_y / lim;
  } else {
    lim = blocks_number_y;
    sm_y = r;
    sm_x = length_x / lim;
  }
  const int sm_x_ = sm_x;
  const int sm_y_ = sm_y;
  double time = 0, start, end;
  time = 0;
  const blitz::Range x_block_first(1, sm_x_), y_block_first(1, sm_y_);
  blitz::Array<blitz::Range, 1> Block_ranges_x(lim), Block_ranges_y(lim);
  for (int v = 1; v <= lim; ++v) {
    Block_ranges_x(v - 1) = x_block_first + (sm_x_ * (v - 1));
    Block_ranges_y(v - 1) = y_block_first + (sm_y_ * (v - 1));
  }
  // array of blocks inner ranges instead
  for (int l = 0; l < repeats; ++l) {
    data(0)(kEz, i, j) = initial;
    data(0)(kHx, i, j) = initial;
    data(0)(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j) = src;
      for (int v1 = 1; v1 <= lim; ++v1) {
       for (int v2 = 1; v2 <= lim; ++v2) {
         UpdateFields2DStencil(Block_ranges_x(v1 - 1), Block_ranges_y(v2 - 1));
       }
    }
      cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
}
  if (time_flag == 1) {
    printf("\n>> Stencil 2D(%d blocks) took %f \n", lim,
           static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    printf("\n>> Stencil 2D (%d blocks) control sum for Ez is ", lim);
    std::cout << sum(abs(data(0)(static_cast<int>(kEz), i, j))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kHx), i, j);
    std::cout << data(0)(static_cast<int>(kHy), i, j);
    std::cout << data(0)(static_cast<int>(kEz), i, j);
  }
  return (time / repeats);
}

void UpdateFields2DRange(blitz::Range temp_x, blitz::Range temp_y) {
  data(1)(kHx, temp_x, temp_y) = data(0)(kChxh, temp_x, temp_y) *
    data(0)(kHx, temp_x, temp_y) - data(0)(kChxe, temp_x, temp_y) *
    (data(0)(kEz, temp_x, temp_y+1) - data(0)(kEz, temp_x, temp_y));

  data(1)(kHy, temp_x, temp_y) = data(0)(kChyh, temp_x, temp_y) *
    data(0)(kHy, temp_x, temp_y) + data(0)(kChye, temp_x, temp_y) *
    (data(0)(kEz, temp_x+1, temp_y) - data(0)(kEz, temp_x, temp_y));

  data(1)(kEz, temp_x, temp_y) = data(0)(kCeze, temp_x, temp_y) *
    data(0)(kEz, temp_x, temp_y) + data(0)(kCezh, temp_x, temp_y) *
    ((data(0)(kHy, temp_x, temp_y) - data(0)(kHy, temp_x-1, temp_y))
     - (data(0)(kHx, temp_x, temp_y) - data(0)(kHx, temp_x, temp_y-1)))
    + data(0)(kSrcEz, temp_x, temp_y);
}
double UsingRange(int length_x, int length_y, int steps) {
  const blitz::Range inner_i(1, length_x);
  const blitz::Range inner_j(1, length_y);
  double time = 0, start, end;
  for (int l = 0; l < repeats; ++l) {
    data(0)(kEz, i, j) = initial;
    data(0)(kHx, i, j) = initial;
    data(0)(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =   exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j) = src;
      UpdateFields2DRange(inner_i, inner_j);
      cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> Range 2D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    std::cout << "\n>> Range 2D control sum for Ez is "
              << sum(abs(data(0)(static_cast<int>(kEz), i, j))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kHx), i, j);
    std::cout << data(0)(static_cast<int>(kHy), i, j);
    std::cout << data(0)(static_cast<int>(kEz), i, j);
  }
  return (time / repeats);
}

double RunAllBlocksForModelRange(int length_x, int length_y, int steps,
                            int r) {
  int blocks_number_y = ceil(length_y / r);
  int blocks_number_x = ceil(length_x / r);
  int lim, sm_x, sm_y;
  if (blocks_number_x > blocks_number_y) {
    lim = blocks_number_x;
    sm_x = r;
    sm_y = length_y / lim;
  } else {
    lim = blocks_number_y;
    sm_y = r;
    sm_x = length_x / lim;
  }
  const int sm_x_ = sm_x;
  const int sm_y_ = sm_y;
  double time = 0, start, end;
  time = 0;
  const blitz::Range x_block_first(1, sm_x_), y_block_first(1, sm_y_);
  blitz::Array<blitz::Range, 1> Block_ranges_x(lim), Block_ranges_y(lim);
  for (int v = 1; v <= lim; ++v) {
    Block_ranges_x(v - 1) = x_block_first + (sm_x_ * (v - 1));
    Block_ranges_y(v - 1) = y_block_first + (sm_y_ * (v - 1));
  }
  for (int l = 0; l < repeats; ++l) {
    data(0)(kEz, i, j) = initial;
    data(0)(kHx, i, j) = initial;
    data(0)(kHy, i, j) = initial;
    start = MPI_Wtime();
    for (int t = 0; t < steps; ++t) {
      const double src =  exp(-static_cast<double>(t * t) / (steps * steps));
      data(0)(static_cast<int>(kSrcEz), i, j) = src;
      for (int v1 = 1; v1 <= lim; ++v1) {
       for (int v2 = 1; v2 <= lim; ++v2) {
         UpdateFields2DRange(Block_ranges_x(v1 - 1), Block_ranges_y(v2 - 1));
       }
    }
      cycleArrays(data(0), data(1));
    }
    end = MPI_Wtime();
    time += end - start;
}
  if (time_flag == 1) {
    printf("\n>> Range 2D(blocks) took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    printf("\n>> Range 2D(%d blocks) control sum for Ez is ", lim);
    std::cout << sum(abs(data(0)(static_cast<int>(kEz), i, j))) << std::endl;
  }
  if (print_flag == 1) {
    std::cout << data(0)(static_cast<int>(kHx), i, j);
    std::cout << data(0)(static_cast<int>(kHy), i, j);
    std::cout << data(0)(static_cast<int>(kEz), i, j);
  }
  return (time / repeats);
}

void PreparePlotDataFile(int length_x, int length_y, int steps) {
  char basename[80] = "blitz_size_", filename[100];
  int max = (length_x < length_y) ? length_y:length_x;
  if (max < 100) {
    snprintf(filename, sizeof(filename), "%s00%dx00%d.plot",
             basename, length_x, length_y);
  } else {
    if (max < 1000) {
      snprintf(filename, sizeof(filename), "%s0%dx0%d.plot",
               basename, length_x, length_y);
    } else {
      snprintf(filename, sizeof(filename), "%s%dx%d.plot",
               basename, length_x, length_y);
    }
  }
  FILE *f;
  f = fopen(filename, "a");
  fprintf(f, "#Model dimensions: %dx%dx%d(steps)\n",
          length_x, length_y, steps);
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

void WritePlotData(int length_x, int length_y, int blocks,
                   double t1, double t2, double t3, double t4, int steps) {
  char basename[80] = "blitz_size_", filename[100];
  int max = (length_x < length_y) ? length_y:length_x;
  if (max < 100) {
    snprintf(filename, sizeof(filename), "%s00%dx00%d.plot",
             basename, length_x, length_y);
  } else {
    if (max < 1000) {
      snprintf(filename, sizeof(filename), "%s0%dx0%d.plot",
               basename, length_x, length_y);
    } else {
      snprintf(filename, sizeof(filename), "%s%dx%d.plot",
               basename, length_x, length_y);
    }
  }
  FILE *f;
  f = fopen(filename, "a");
  if (block_flag == 0) {
    int total = length_x * length_y * steps;
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
  int total = length_x * length_y * steps;
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
  int C = max_steps * max_size * max_size;
  int counter1 = 0, counter2;
  time_stencil.resize(1 + GetPowerOfTwo(max_size / min_size));
  time_stencil_blocks.resize((1 + GetPowerOfTwo(max_size / min_size)),
                           (GetPowerOfTwo(sm_max / sm_min)));
  time_range.resize(1 + GetPowerOfTwo(max_size / min_size));
  time_range_blocks.resize((1 + GetPowerOfTwo(max_size / min_size)),
                           (GetPowerOfTwo(sm_max / sm_min)));
  //  time_range_blocks = 0.0;
  blitz::Array <int, 1> blocks(GetPowerOfTwo(sm_max / sm_min));

  for (int size = min_size; size <= max_size; size *=2) {
    length_x = size;
    length_y = size;
    steps = C / (length_x * length_y);
    PreparePlotDataFile(length_x, length_y, steps);
    for (int counter = 0; counter < depth; ++counter) {
      data(counter).resize(components, length_x + 2, length_y + 2);
      data(counter) = 0;
    }
    double coefficient = 0.1;
    for (int w = 0; w < depth; w++) {
      data(w)(kCeze, i, j) = 1;
      data(w)(kCezh, i, j) = coefficient;
      data(w)(kChxe, i, j) = coefficient;
      data(w)(kChxh, i, j) = coefficient;
      data(w)(kChye, i, j) = coefficient;
      data(w)(kChyh, i, j) = coefficient;
    }
      printf("\n>> __Processing task: 2D %dx%dx%d(steps)\n",
             length_x, length_y, steps);
    time_stencil(counter1) = UsingStencils(length_x, length_y, steps);
    time_range(counter1) = UsingRange(length_x, length_y, steps);

    int max = (length_x <= length_y) ? length_y : length_x;
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
            RunAllBlocksForModelStencil(length_x, length_y, steps, r);
          time_range_blocks(counter1, counter2) =
          RunAllBlocksForModelRange(length_x, length_y, steps, r);
          }
        }
        counter2 += 1;
      }
    }  // end of if (sm_min < max)
    for (int cnt = 0; cnt < (counter2 - 1); ++cnt) {
      WritePlotData(length_x, length_y, blocks(cnt),
                    time_stencil(counter1), time_stencil_blocks(counter1, cnt),
                    time_range(counter1), time_range_blocks(counter1, cnt), steps);
    }
    counter1 += 1;
      printf("\n>> __Done processing task: 2D %dx%dx%d(steps)\n",
             length_x, length_y, steps);
  }  // end of for(int size = min_size; size <= max_size; size *= 2)
}
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int total_processes_number, process_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  // Output parameters **********************************************************
  sum_flag = 1;
  print_flag = 0;
  time_flag = 1;
  block_flag = blk;
  int size_check_flag = 1;
  // ***************************************************************************
  // Parameters to set *********************************************************
  // set number of timesteps and number of time measure procedure repeats
  max_steps = 10;
  repeats = 2;
  // set min ( >= 16 ) and max size of model
  min_size = 32;
  max_size = 2048;
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

