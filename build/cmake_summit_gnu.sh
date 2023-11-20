#!/bin/bash

source ${MODULESHOME}/init/bash
module purge
module load DefApps gcc/9.3.0 cuda parallel-netcdf cmake

unset CUDAFLAGS
unset CXXFLAGS

./cmake_clean.sh

cmake -DCMAKE_CXX_COMPILER=mpic++                   \
      -DCMAKE_C_COMPILER=mpicc                      \
      -DCMAKE_Fortran_COMPILER=mpif90               \
      -DYAKL_ARCH="CUDA"                            \
      -DYAKL_CUDA_FLAGS="-DHAVE_MPI -O3 --use_fast_math -arch sm_70 -ccbin mpic++ -I${OLCF_PARALLEL_NETCDF_ROOT}/include" \
      -DLDFLAGS="-L${OLCF_PARALLEL_NETCDF_ROOT}/lib -lpnetcdf"  \
      ..

