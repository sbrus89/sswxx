#!/bin/bash

source ${MODULESHOME}/init/bash
module purge
module load DefApps gcc/9.3.0 cuda parallel-netcdf cmake

unset CUDAFLAGS
unset CXXFLAGS
