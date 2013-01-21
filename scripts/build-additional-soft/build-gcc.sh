#!/bin/bash
cores=4  # should be common for most users going to build gcc...
path_build_gcc=$PWD/gcc-4.7
path_sources=$path_build_gcc/sources
path_build=$path_build_gcc/build
path_output=$path_build_gcc/output
ppl=ppl-0.11.2
cloog=cloog-0.16.1
mpc=mpc-1.0.1
gmp=gmp-5.0.5
mpfr=mpfr-3.1.1
gcc=gcc-4.7.2
file_ppl=$ppl.tar.bz2
file_cloog=$cloog.tar.gz
file_mpc=$mpc.tar.gz
file_gmp=$gmp.tar.bz2
file_mpfr=$mpfr.tar.bz2
file_gcc=$gcc.tar.bz2
echo Prepare to build GCC 4.7 inside $path_build_gcc
if [[ -a $path_build_gcc ]]; then
    echo Found build folder.
else
    echo Creating build folder...
    mkdir $path_build_gcc
fi
###########################################
echo Download sources...
if [[ -a $path_sources ]]; then
    echo Found sources folder.
else
    echo Creating sources folder...
    mkdir $path_sources
fi
cd $path_sources
# if [[ -a $file_ppl ]]; then
#     echo Found PPL $file_ppl
# else
#     echo Downloading PPL
#     # wget http://194.85.163.135/$file_ppl
#     wget http://pkgs.fedoraproject.org/repo/pkgs/ppl/ppl-0.11.2.tar.bz2/c24429e6c3bc97d45976a63f40f489a1/ppl-0.11.2.tar.bz2
# fi
# if [[ -a $file_cloog ]]; then
#     echo Found CLooG $file_cloog
# else
#     echo Downloading CLooG
#     wget http://www.bastoul.net/cloog/pages/download/$file_cloog
# fi
# if [[ -a $file_mpc ]]; then
#     echo Found MPC $file_mpc
# else
#     echo Downloading MPC
#     wget http://www.multiprecision.org/mpc/download/$file_mpc
# fi
# if [[ -a $file_gmp ]]; then
#     echo Found GMP $file_gmp
# else
#     echo Downloading GMP
#     wget ftp://ftp.gmplib.org/pub/${file_gmp/.tar.bz2/}/$file_gmp
# fi
# if [[ -a $file_mpfr ]]; then
#     echo Found MPFR $file_mpfr
# else
#     echo Downloading MPFR
#     wget http://mpfr.loria.fr/mpfr-current/$file_mpfr
# fi
if [[ -a $file_gcc ]]; then
    echo Found GCC $file_gcc
else
    echo Downloading GCC
    wget http://www.netgull.com/gcc/releases/${file_gcc/.tar.bz2/}/$file_gcc
fi
###########################################
echo Unpack sources...
# for file in *.tar.gz; do
#     echo Unpacking $file ...
#     tar xzf $file
# done
for file in *.tar.bz2; do
    echo Unpacking $file ...
    tar xjf $file
done

# ###########################################
echo Unpack sources...
if [[ -a $path_build ]]; then
    echo Found build folder.
else
    echo Creating build folder...
    mkdir $path_build
fi
if [[ -a $path_output ]]; then
    echo Found output folder.
else
    echo Creating output folder...
    mkdir $path_output
fi
# ###########################################
# echo Build...
# echo Build GMP: 
# if [[ ! -a $path_build/$gmp ]]; then
#     mkdir $path_build/$gmp
# fi
# cd $path_build/$gmp
# $path_sources/$gmp/configure --prefix=$path_output/
# make && make check && make install

# echo Build MPFR: 
# if [[ ! -a $path_build/$mpfr ]]; then
#     mkdir $path_build/$mpfr
# fi
# cd $path_build/$mpfr
# $path_sources/$mpfr/configure \
#     --prefix=$path_output/ \
#     --with-gmp=$path_output/
# make && make install

# echo Build MPC: 
# if [[ ! -a $path_build/$mpc ]]; then
#     mkdir $path_build/$mpc
# fi
# cd $path_build/$mpc
# $path_sources/$mpc/configure \
#     --prefix=$path_output/ \
#     --with-gmp=$path_output/ \
#     --with-mpfr=$path_output/
# make && make install

# echo Build PPL: 
# if [[ ! -a $path_build/$ppl ]]; then
#     mkdir $path_build/$ppl
# fi
# cd $path_build/$ppl
# $path_sources/$ppl/configure \
#     --prefix=$path_output/ \
#     --with-gmp=$path_output/
# make && make install

# echo Build CLooG: 
# if [[ ! -a $path_build/$cloog ]]; then
#     mkdir $path_build/$cloog
# fi
# cd $path_build/$cloog
# $path_sources/$cloog/configure \
#     --prefix=$path_output/ \
#     --with-gmp=$path_output/
# make && make install

echo Build GCC: 
if [[ ! -a $path_build/$gcc ]]; then
    mkdir $path_build/$gcc
fi
cd $path_sources/${file_gcc/.tar.bz2/}
./contrib/download_prerequisites
cd $path_build/$gcc
make distclean
# $path_sources/$gcc/configure \
#     --prefix=$path_output/    \
#     --enable-languages=c,c++,fortran 

$path_sources/$gcc/configure \
    --prefix=$path_output/    \
    --enable-languages=c,c++,fortran,objc,obj-c++ \
    --program-suffix=-4.7 \
    --enable-shared \
    --enable-multiarch  \
    --enable-linker-build-id \
#    --with-system-zlib \
    --without-included-gettext \
    --enable-threads=posix \
    --enable-nls  \
    --enable-clocale=gnu \
    --enable-libstdcxx-debug \
    --enable-objc-gc \
    --with-arch-32=i586 \
    --with-tune=generic \
    --enable-checking=release \
    --build=x86_64-linux-gnu \
    --host=x86_64-linux-gnu \
    --target=x86_64-linux-gnu

# $path_sources/$gcc/configure \
#     --prefix=$path_output/    \
#     --with-gmp=$path_output/  \
#     --with-mpfr=$path_output/ \
#     --with-mpc=$path_output/  \
#     --enable-cloog-backend=isl       \
#     --enable-languages=c,c++,fortran \
#     --enable-lto \
#     --with-ppl=$path_output/  \
#     --with-cloog=$path_output/

export LD_LIBRARY_PATH=$path_output/lib/
make && make install
