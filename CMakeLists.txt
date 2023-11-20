cmake_minimum_required(VERSION 3.8)
project(ssw CXX)

if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()


############################################################
## Append CXXFLAGS
############################################################
SET(CMAKE_CXX_FLAGS "${CXXFLAGS}")


############################################################
## Add YAKL 
############################################################
add_subdirectory(YAKL)
if (${CMAKE_VERSION} VERSION_GREATER "3.18.0")
  set_property(TARGET yakl PROPERTY CUDA_ARCHITECTURES OFF)
endif()
include_directories(YAKL)


############################################################
## YAKL parallel_for Version
############################################################
add_executable(ssw_model driver.cpp)
set_target_properties(ssw_model PROPERTIES COMPILE_FLAGS "${EXE_DEFS}")

target_link_libraries(ssw_model "${LDFLAGS}")


include(YAKL/yakl_utils.cmake)
yakl_process_target(ssw_model)

if ("${YAKL_ARCH}" STREQUAL "CUDA")
  set_target_properties(ssw_model PROPERTIES LINKER_LANGUAGE CXX)
  if (${CMAKE_VERSION} VERSION_GREATER "3.18.0")
    set_target_properties(ssw_model PROPERTIES CUDA_ARCHITECTURES OFF)
  endif()
endif()



