#!/bin/bash

source ${MODULESHOME}/init/bash
module purge
#module load DefApps gcc/9.3.0 cuda parallel-netcdf cmake
module load DefApps gcc/9.3.0 cuda cmake netcdf-c/4.8.1 parallel-netcdf

unset CUDAFLAGS
unset CXXFLAGS

./cmake_clean.sh

cmake -DCMAKE_CXX_COMPILER=mpic++                   \
      -DCMAKE_C_COMPILER=mpicc                      \
      -DCMAKE_Fortran_COMPILER=mpif90               \
      -DYAKL_ARCH=""                            \
      -DCXXFLAGS="-O0 -g -DARRAY_DEBUG -DYAKL_DEBUG"                            \
      -DYAKL_CUDA_FLAGS="-DHAVE_MPI -DARRAY_DEBUG -O0 -g  -arch sm_70 -ccbin mpic++ -I${OLCF_PARALLEL_NETCDF_ROOT}/include -I${OLCF_NETCDF_C_ROOT}/include" \
      -DLDFLAGS="-L${OLCF_NETCDF_C_ROOT}/lib -lnetcdf -L${OLCF_PARALLEL_NETCDF_ROOT}/lib -lpnetcdf"  \
      ..

