#!/bin/bash
# Running Meep
file="task2D.cc"
echo -e "\n >> Meep 2D\n"
g++ -malign-double $file -DOnzaTestSize=2048 -DOnzaTestSteps=100 -o $file.bin /home/dmmrkovich/meep-inst/meep-1.1.1/src/.libs/libmeep.a -lhdf5 -lz -lgsl -lharminv -llapack -latlas -lfftw3 -lm
#-lcblas
./$file.bin
file="task3D.cc"
echo -e "\n >> Meep 3D\n"
g++ -malign-double $file -DOnzaTestSize=128 -DOnzaTestSteps=100 -o $file.bin /home/dmmrkovich/meep-inst/meep-1.1.1/src/.libs/libmeep.a -lhdf5 -lz -lgsl -lharminv -llapack -latlas -lfftw3 -lm
#-lcblas
./$file.bin
# End of Running Meep
# Running StdC
file="2D_test_c_meep.c"
echo -e "\n >> StdC 2D\n"
mpicc -O3 -Dsz=2048 -Dst=100 -fstrict-aliasing -std=c99 -lm $file -o $file.bin
./$file.bin
file="3D_test_c_meep.c"
echo -e "\n >> StdC 3D\n"
mpicc -O3 -Dsz=128 -Dst=100 -fstrict-aliasing -std=c99 -lm $file -o $file.bin
./$file.bin
# End of Running StdC
rm -r *.bin
