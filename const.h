
#pragma once

#include "YAKL.h"
#include <cmath>

using yakl::SArray;
using yakl::c::SimpleBounds;

#ifdef SINGLE_PREC
  typedef float  real;
  auto mpi_type = MPI_FLOAT;
#else
  typedef double real;
  auto mpi_type = MPI_DOUBLE;
#endif

using yakl::c::Bounds;
using yakl::c::parallel_for;
using yakl::SArray;

template<class T> inline T min( T val1 , T val2 ) {
  return val1 < val2 ? val1 : val2 ;
}

template<class T> inline T abs( T val ) {
  return val > 0 ? val : -val;
}

#ifdef SIMD_LEN
  unsigned int constexpr simd_len = SIMD_LEN;
#else
  unsigned int constexpr simd_len = 4;
#endif

using yakl::simd::Pack;
using yakl::simd::PackIterConfig;
using yakl::simd::iterate_over_pack;


typedef yakl::Array<int   ,1,yakl::memDevice> int1d;
typedef yakl::Array<int   ,2,yakl::memDevice> int2d;
typedef yakl::Array<int   ,3,yakl::memDevice> int3d;
typedef yakl::Array<real  ,1,yakl::memDevice> real1d;
typedef yakl::Array<real  ,2,yakl::memDevice> real2d;
typedef yakl::Array<real  ,3,yakl::memDevice> real3d;
typedef yakl::Array<double,1,yakl::memDevice> doub1d;
typedef yakl::Array<double,2,yakl::memDevice> doub2d;
typedef yakl::Array<double,3,yakl::memDevice> doub3d;

typedef yakl::Array<int    const,1,yakl::memDevice> intConst1d;
typedef yakl::Array<int    const,2,yakl::memDevice> intConst2d;
typedef yakl::Array<int    const,3,yakl::memDevice> intConst3d;
typedef yakl::Array<real   const,1,yakl::memDevice> realConst1d;
typedef yakl::Array<real   const,2,yakl::memDevice> realConst2d;
typedef yakl::Array<real   const,3,yakl::memDevice> realConst3d;
typedef yakl::Array<double const,1,yakl::memDevice> doubConst1d;
typedef yakl::Array<double const,2,yakl::memDevice> doubConst2d;
typedef yakl::Array<double const,3,yakl::memDevice> doubConst3d;

typedef yakl::Array<int   ,1,yakl::memHost> int1dHost;
typedef yakl::Array<int   ,2,yakl::memHost> int2dHost;
typedef yakl::Array<int   ,3,yakl::memHost> int3dHost;
typedef yakl::Array<real  ,1,yakl::memHost> real1dHost;
typedef yakl::Array<real  ,2,yakl::memHost> real2dHost;
typedef yakl::Array<real  ,3,yakl::memHost> real3dHost;
typedef yakl::Array<double,1,yakl::memHost> doub1dHost;
typedef yakl::Array<double,2,yakl::memHost> doub2dHost;
typedef yakl::Array<double,3,yakl::memHost> doub3dHost;
