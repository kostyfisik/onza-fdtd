/// @file   3D_test_std.cc
/// @author Markovich Dmitry <dmmrkovich at gmail (.) com>
/// @copyright 2012 Markovich Dmitry
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <cstdlib>
#define SM (CLS / sizeof(double))

const int depth = 2;
int repeats, min_size, max_size, max_steps;
int sum_flag, print_flag, time_flag, size_check_flag;
double* time_std_c;
double** time_std_c_blocks;
const double initial = 0.0, coefficient = 0.1;
enum DataComponents {kEx = 0, kCexe, kCexh,
                     kEy, kCeye, kCeyh,
                     kEz, kCeze, kCezh,
                     kHx, kChxh, kChxe,
                     kHy, kChyh, kChye,
                     kHz, kChzh, kChze,
                     kSrcEz };
const int components = kSrcEz + 1;
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

class FDTD {
  public:
  int length_x_, length_y_, length_z_, steps_;
  double *****data_;
  FDTD(int length_x, int length_y, int length_z, int steps, int repeats);
  ~FDTD();
  void UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2);
  void InitializeFields();
  void SetCoefficients();
  void CycleArrays();
  double EzControlSum();
  void PrintField();
  double UsingStdC();
  double UsingStdCBlocks(int block_length);
};

FDTD::FDTD(int length_x, int length_y, int length_z, int steps, int repeats) {
  length_x_ = length_x;
  length_y_ = length_y;
  length_z_ = length_z;
  steps_ = steps;

  data_ = new double**** [depth];

  for (int i1 = 0; i1 < depth; ++i1) {
    data_[i1] = new double***[components];
  }
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i2 = 0; i2 < components; ++i2) {
      data_[i1][i2] = new double**[length_x_ + 2];
    }
  }
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i2 = 0; i2 < components; ++i2) {
      for (int i3 = 0; i3 <= length_x + 1; ++i3) {
      data_[i1][i2][i3] = new double*[length_y_ + 2];
      }
    }
  }
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i2 = 0; i2 < components; ++i2) {
      for (int i3 = 0; i3 <= length_x + 1; ++i3) {
	for (int i4 = 0; i4 <= length_y + 1; ++i4) {
	  data_[i1][i2][i3][i4] = new double[length_z_ + 2];
	}
      }
    }
  }

  printf("\n>> __Now processing task 3D %dx%dx%dx%d(steps)\n", length_x_,
	 length_y_, length_z_, steps_);
}
double ABS(double number) {
  if (number >= 0) {
    return number;
  } else {
    return -1 * number;
  }
}
FDTD::~FDTD() {
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i2 = 0; i2 < components; ++i2) {
      for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
	for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
	  delete[]data_[i1][i2][i3][i4];
	}
      }
    }
  }
  for (int i1 = 0; i1 < depth; ++i1) {
    delete[]data_[i1];
  }
  delete[] data_;
printf("\n>> __Done processing task 3D %dx%dx%d(steps)\n", length_x_, length_y_,
       length_z_, steps_);
}
void FDTD::UpdateFields3DStdC(int x1, int x2, int y1, int y2, int z1, int z2) {
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kHx][i3][i4][i5] = data_[0][kChxh][i3][i4][i5] *
	  data_[0][kHx][i3][i4][i5] + data_[0][kChxe][i3][i4][i5] *
	  ((data_[0][kEy][i3][i4][i5 + 1] - data_[0][kEy][i3][i4][i5]) -
	   (data_[0][kEz][i3][i4 + 1][i5] - data_[0][kEz][i3][i4][i5]));
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kHy][i3][i4][i5] = data_[0][kChyh][i3][i4][i5] *
	  data_[0][kHy][i3][i4][i5] + data_[0][kChye][i3][i4][i5] *
	  ((data_[0][kEz][i3 + 1][i4][i5] - data_[0][kEz][i3][i4][i5])
	   -(data_[0][kEx][i3][i4][i5 + 1] - data_[0][kEx][i3][i4][i5]));
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kHz][i3][i4][i5] = data_[0][kChzh][i3][i4][i5] *
	  data_[0][kHz][i3][i4][i5] + data_[0][kChze][i3][i4][i5] *
	  ((data_[0][kEx][i3][i4 + 1][i5] - data_[0][kEx][i3][i4][i5])
	   -(data_[0][kEy][i3 + 1][i4][i5] - data_[0][kEy][i3][i4][i5]));
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kEx][i3][i4][i5] = data_[0][kCexe][i3][i4][i5] *
	  data_[0][kEx][i3][i4][i5] + data_[0][kCexh][i3][i4][i5] *
	  ((data_[0][kHz][i3][i4][i5] - data_[0][kHz][i3][i4 - 1][i5])
	   - (data_[0][kHy][i3][i4][i5] - data_[0][kHy][i3][i4][i5 - 1]));
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kEy][i3][i4][i5] = data_[0][kCeye][i3][i4][i5] *
	  data_[0][kEy][i3][i4][i5] + data_[0][kCeyh][i3][i4][i5] *
	  ((data_[0][kHx][i3][i4][i5] - data_[0][kHx][i3][i4][i5 - 1])
	   - (data_[0][kHz][i3][i4][i5] - data_[0][kHz][i3 - 1][i4][i5]));
      }
    }
  }
  for (int i3 = x1; i3 <= x2; ++i3) {
    for (int i4 = y1; i4 <= y2; ++i4) {
      for (int i5 = z1; i5 <= z2; ++i5) {
	data_[1][kEz][i3][i4][i5] = data_[0][kCeze][i3][i4][i5] *
	  data_[0][kEz][i3][i4][i5] + data_[0][kCezh][i3][i4][i5] *
	  ((data_[0][kHy][i3][i4][i5] - data_[0][kHy][i3 - 1][i4][i5])
	   - (data_[0][kHx][i3][i4][i5] - data_[0][kHx][i3][i4 - 1][i5]))
	  + data_[0][kSrcEz][i3][i4][i5];
      }
    }
  }
}
void FDTD::InitializeFields() {
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
      for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
	for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	  data_[i1][kHx][i3][i4][i5] = initial;
	  data_[i1][kHy][i3][i4][i5] = initial;
	  data_[i1][kHz][i3][i4][i5] = initial;
	  data_[i1][kEx][i3][i4][i5] = initial;
	  data_[i1][kEy][i3][i4][i5] = initial;
	  data_[i1][kEz][i3][i4][i5] = initial;
	  data_[i1][kSrcEz][i3][i4][i5] = initial;
	}
      }
    }
  }
}
void FDTD::SetCoefficients() {
  for (int i1 = 0; i1 < depth; ++i1) {
    for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
      for (int i4 = 0; i4 <= length_y_+ 1; ++i4) {
	for (int i5 = 0; i5 <= length_z_+ 1; ++i5) {
	  data_[i1][kChxh][i3][i4][i5] = coefficient;
	  data_[i1][kChxe][i3][i4][i5] = coefficient;
	  data_[i1][kChyh][i3][i4][i5] = coefficient;
	  data_[i1][kChye][i3][i4][i5] = coefficient;
	  data_[i1][kChzh][i3][i4][i5] = coefficient;
	  data_[i1][kChze][i3][i4][i5] = coefficient;
	  data_[i1][kCexe][i3][i4][i5] = coefficient;
	  data_[i1][kCexh][i3][i4][i5] = coefficient;
	  data_[i1][kCeye][i3][i4][i5] = coefficient;
	  data_[i1][kCeyh][i3][i4][i5] = coefficient;
	  data_[i1][kCeze][i3][i4][i5] = 1;
	  data_[i1][kCezh][i3][i4][i5] = coefficient;
	}
      }
    }
  }
}

void FDTD::CycleArrays() {
  for (int i2 = 0; i2 < components; ++i2) {
    for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
      for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
	for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	  data_[0][i2][i3][i4][i5] = data_[1][i2][i3][i4][i5];
	}
      }
    }
  }
}

double FDTD::EzControlSum() {
  double sum = 0;
  for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	sum += ABS(data_[0][kEz][i3][i4][i5]);
      }
    }
  }
  return sum;
}
void FDTD::PrintField() {
  for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
    for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
      for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	printf("%f\t", data_[0][kEz][i3][i4][i5]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}
double FDTD::UsingStdC() {
  double time = 0, start, end;
  SetCoefficients();
  for (int o = 0; o < repeats; ++o) {
    InitializeFields();
    start = MPI_Wtime();
    for (int t = 0; t < steps_; ++t) {
      for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
        for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
	  for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	    data_[0][kSrcEz][i3][i4][i5] =
	      exp(-static_cast<double>(t * t) / (steps_ * steps_));
	  }
	}
      }
      UpdateFields3DStdC(1, length_x_, 1, length_y_, 1, length_z_);
      CycleArrays();
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> StdC 3D took %f \n", static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    printf("\n>> StdC 3D Ez control sum is %e\n",  EzControlSum());
  }
  if (print_flag == 1) {
    PrintField();
  }
  return (time / repeats);
}

double FDTD::UsingStdCBlocks(int block_length) {
  int max = (length_x_ <= length_y_) ? length_y_:length_x_;
  max = (max <= length_z_) ? length_z_:max;
  int lim = ceil(max / block_length);
  const int sm_x = length_x_ / lim, sm_y = length_y_ / lim,
    sm_z = length_z_ / lim;
  double time = 0, start, end;
  time = 0;
  int Block_ranges_x1[lim], Block_ranges_x2[lim],
    Block_ranges_y1[lim], Block_ranges_y2[lim],
    Block_ranges_z1[lim], Block_ranges_z2[lim];
  // printf("\n>> Using %d blocks, sm_x = %d, sm_y = %d\n", lim, sm_x, sm_y);
  for (int v = 1; v <= lim; ++v) {
    Block_ranges_x1[v] = 1 + sm_x * (v - 1);
    Block_ranges_x2[v] = sm_x + sm_x * (v - 1);
    Block_ranges_y1[v] = 1 + sm_y * (v - 1);
    Block_ranges_y2[v] = sm_y + sm_y * (v - 1);
    Block_ranges_z1[v] = 1 + sm_z * (v - 1);
    Block_ranges_z2[v] = sm_z + sm_z * (v - 1);
  }
  for (int l = 0; l < repeats; ++l) {
    InitializeFields();
    start = MPI_Wtime();
    for (int t = 0; t < steps_; ++t) {
      for (int i3 = 0; i3 <= length_x_ + 1; ++i3) {
        for (int i4 = 0; i4 <= length_y_ + 1; ++i4) {
	  for (int i5 = 0; i5 <= length_z_ + 1; ++i5) {
	    data_[0][kSrcEz][i3][i4][i5] =
	      exp(-static_cast<double>(t * t) / (steps_ * steps_));
	  }
        }
      }
      for (int v1 = 1; v1 <= lim; ++v1) {
        for (int v2 = 1; v2 <= lim; ++v2) {
	  for (int v3 = 1; v3 <= lim; ++v3) {
	    UpdateFields3DStdC(Block_ranges_x1[v1],
			       Block_ranges_x2[v1],
			       Block_ranges_y1[v2],
			       Block_ranges_y2[v2],
			       Block_ranges_z1[v3],
			       Block_ranges_z2[v3]);
	  }
	}
      }
      CycleArrays();
    }
    end = MPI_Wtime();
    time += end - start;
  }
  if (time_flag == 1) {
    printf("\n>> StdC 3D (%d blocks) took %f \n", lim,
           static_cast<float>(time / repeats));
  }
  if (sum_flag == 1) {
    printf("\n>> StdC 3D (%d blocks) Ez control sum is %e\n",  lim, EzControlSum());
  }
  if (print_flag == 1) {
    PrintField();
  }
  return (time / repeats);
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

void PreparePlotDataFile(int length_x, int length_y, int length_z, int steps) {
  char basename[80] = "std_c_results_size_", filename[100];
  int max = (length_x < length_y) ? length_y:length_x;
  max = (max <= length_z) ? length_z:max;
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
  fprintf(f, "#Model dimensions: %dx%dx%dx%d(steps)\n", length_x, length_y, 
	  length_z, steps);
  fprintf(f, "#Blocks' length\t\tStdC\t\tStdC(blocks)\t\tBest");
  fprintf(f, "\t\tt(std_c)\t\tt(std_c_blcks)\n");
  fclose(f);
}

void WritePlotData(int length_x, int length_y, int length_z, int blocks,
                   double t1, double t2, int steps) {
  char basename[80] = "std_c_results_size_", filename[100];
  int max = (length_x <= length_y) ? length_y:length_x;
  max = (length_z <= max) ? max:length_z;
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
  fprintf(f, "%d\t\t\t%.2f\t\t%.2f\t", blocks, t1, t2);
  double min = t1;
  min = (t2 < min) ? t2:min;
  if (min == t1) {
    fprintf(f, "\tStdC");
  }
  if (min == t2) {
    fprintf(f, "\tStdC(blocks)(+%.1f%)", (100 * (t1-t2) / t1));
  }
  int total = length_x * length_y * length_z * steps;
  fprintf(f, "\t%.2e\t\t%.2e\n", (t1 / total), (t2 / total));
  fclose(f);
}

void RunCalculation() {
  int sm_min = CheckCLS(), sm_max = max_size;
  int length_x, length_y, length_z, steps;
  int C = max_size * max_size * max_size * max_steps;
  const int kDim1 = 1 + GetPowerOfTwo(max_size / min_size),
    kDim2 = GetPowerOfTwo(sm_max / sm_min);
  time_std_c = new double[kDim1];
  time_std_c_blocks = new double*[kDim1];
  for (int i = 0; i < (1 + GetPowerOfTwo(max_size / min_size)); ++i) {
    time_std_c_blocks[i] = new double[kDim2];
  }
  int blocks[kDim2];
  int counter1 = 0;
  for (int size = min_size; size <= max_size; size *= 2) {
    length_x = size;
    length_y = size;
    length_z = size;
    steps = C / (length_x * length_y * length_z);
    PreparePlotDataFile(length_x, length_y, length_z, steps);
    FDTD task(length_x, length_y, length_z, steps, repeats);
    time_std_c[counter1] = task.UsingStdC();
    int max = (length_x <= length_y) ? length_y : length_x;
    max = (max <= length_z) ? length_z : max;
    int counter2 = 0;
    if (sm_min < max) {
      for (int block_length = sm_min; block_length <= max; block_length *= 2) {
        if (block_length < max) {
          blocks[counter2] = block_length;
          time_std_c_blocks[counter1][counter2] =
            task.UsingStdCBlocks(block_length);
        }
        counter2 += 1;
      }
      for (int i2 = 0; i2 < counter2 - 1; ++i2) {
	WritePlotData(length_x, length_y, length_z, blocks[i2],
		      time_std_c[counter1], time_std_c_blocks[counter1][i2],
		      steps);
      }
    }
    counter1 += 1;
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int total_processes_number, process_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  // Flags *********************************************************
  sum_flag = 1;
  print_flag = 0;
  time_flag = 1;
  size_check_flag = 0;
  // ***************************************************************************
  // Parameters to set *********************************************************
  // set number of timesteps and number of time measure procedure repeats
  max_steps = 100;
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
    RunCalculation();
    MPI_Finalize();
    return 0;
}

