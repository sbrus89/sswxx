
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

real const gravity = 9.81;

using yakl::simd::Pack;
using yakl::simd::PackIterConfig;
using yakl::simd::iterate_over_pack;


using int1d = yakl::Array<int   ,1,yakl::memDevice>;
using int2d = yakl::Array<int   ,2,yakl::memDevice>;
using int3d = yakl::Array<int   ,3,yakl::memDevice>;
using real1d = yakl::Array<real  ,1,yakl::memDevice>;
using real2d = yakl::Array<real  ,2,yakl::memDevice>;
using real3d = yakl::Array<real  ,3,yakl::memDevice>;
using doub1d = yakl::Array<double,1,yakl::memDevice>;
using doub2d = yakl::Array<double,2,yakl::memDevice>;
using doub3d = yakl::Array<double,3,yakl::memDevice>;

using intConst1d = yakl::Array<int    const,1,yakl::memDevice>;
using intConst2d = yakl::Array<int    const,2,yakl::memDevice>;
using intConst3d = yakl::Array<int    const,3,yakl::memDevice>;
using realConst1d = yakl::Array<real   const,1,yakl::memDevice>;
using realConst2d = yakl::Array<real   const,2,yakl::memDevice>;
using realConst3d = yakl::Array<real   const,3,yakl::memDevice>;
using doubConst1d = yakl::Array<double const,1,yakl::memDevice>;
using doubConst2d = yakl::Array<double const,2,yakl::memDevice>;
using doubConst3d = yakl::Array<double const,3,yakl::memDevice>;

using int1dHost = yakl::Array<int   ,1,yakl::memHost>;
using int2dHost = yakl::Array<int   ,2,yakl::memHost>;
using int3dHost = yakl::Array<int   ,3,yakl::memHost>;
using real1dHost = yakl::Array<real  ,1,yakl::memHost>;
using real2dHost = yakl::Array<real  ,2,yakl::memHost>;
using real3dHost = yakl::Array<real  ,3,yakl::memHost>;
using doub1dHost = yakl::Array<double,1,yakl::memHost>;
using doub2dHost = yakl::Array<double,2,yakl::memHost>;
using doub3dHost = yakl::Array<double,3,yakl::memHost>;
