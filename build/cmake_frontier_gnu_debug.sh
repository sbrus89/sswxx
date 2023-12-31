#!/bin/bash

source ${MODULESHOME}/init/bash
module purge
#module load PrgEnv-amd cray-parallel-netcdf cmake craype-accel-amd-gfx90a
module load PrgEnv-gnu/8.4.0 cray-parallel-netcdf cmake
#module load PrgEnv-gnu/8.4.0 cmake/3.23.2 cray-hdf5/1.12.2.7 cray-netcdf/4.9.0.7

unset CXX
unset CC
unset FC
unset CUDAFLAGS
unset CXXFLAGS

export MPICH_GPU_SUPPORT_ENABLED=0

./cmake_clean.sh

cmake -DCMAKE_CXX_COMPILER=mpic++   \
      -DCMAKE_C_COMPILER=mpicc                 \
      -DCMAKE_Fortran_COMPILER=mpif90               \
      -DYAKL_ARCH=""                            \
      -DCXXFLAGS="-O0 -g -DARRAY_DEBUG -DYAKL_DEBUG"                            \
      -DYAKL_CXX_FLAGS="-DHAVE_MPI -DARRAY_DEBUG -O0 -g -I${PNETCDF_DIR}/include" \
      -DLDFLAGS="-L${PNETCDF_DIR}/lib -lpnetcdf"  \
      ..

#cmake -DCMAKE_CXX_COMPILER=mpic++   \
#      -DCMAKE_C_COMPILER=mpicc                 \
#      -DCMAKE_Fortran_COMPILER=mpif90               \
#      -DYAKL_ARCH=""                            \
#      -DCXXFLAGS="-O0 -g -DARRAY_DEBUG -DYAKL_DEBUG"                            \
#      -DYAKL_CXX_FLAGS="-DHAVE_MPI -DARRAY_DEBUG -O0 -g -I${NETCDF_DIR}/include" \
#      -DLDFLAGS="-L${NETCDF_DIR}/lib -lnetcdf"  \
#      ..
