#include </home/dmmrkovich/meep-inst/meep-1.1.1/src/meep.hpp>
using namespace meep;

double eps(const vec &p) {
  //if (p.x() < 2 && p.y() < 3)
  //  return 12.0;
  return 1.0;
}

complex<double> one(const vec &p) {
 return 1.0;
}

const int size = OnzaTestSize;
const int steps = OnzaTestSteps;
const double epsilon = 1.0;

int main(int argc, char **argv) {
  initialize mpi(argc, argv); // do this even for non-MPI Meep
 const double resolution = 1; // pixels per distance
  grid_volume v = vol2d(size,size,resolution); // size x size 2d cell
  structure s(v, eps/*, pml(1.0)*/);
  fields f(&s);
  f.use_real_fields();
  
//  f.output_hdf5(Dielectric, v.surroundings());
  
  const double freq = 0.15, fwidth = 0.1;
  const gaussian_src_time src(freq, fwidth);
  const volume src_plane(vec(0,0),vec(size,size));
  f.add_volume_source(Ez,src,src_plane,one,1.0);
// while (f.time() < steps/*f.last_source_time()*/) {
//    f.step();
//  }
   for (int i = 0; i < steps; ++i) {
	f.step();
   }
  
//  f.output_hdf5(Ez, v.surroundings());
  
  return 0;
} 
