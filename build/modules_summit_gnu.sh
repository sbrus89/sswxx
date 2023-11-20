#!/bin/bash

source ${MODULESHOME}/init/bash
module purge
module load DefApps gcc/9.3.0 cuda parallel-netcdf cmake netcdf-c/4.8.1

unset CUDAFLAGS
unset CXXFLAGS
