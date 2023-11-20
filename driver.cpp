#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "const.h"
#include "pnetcdf.h"
//#include "netcdf.h"
#include "io.hpp"
#include "mesh.hpp"
#include "state.hpp"
#include "timestep.hpp"
#include "decomp.hpp"
#include <ctime>
#include <chrono>


int main (int argc, char **argv) {
  real dt, t;
  real end_time;
  IO io;

  end_time = 10.0*3600.0;
  dt = 10.0*60.0;
  t = 0.0;

  MPI_Init(&argc, &argv);
  std::cout << "Initializing yakl" << std::endl;
  yakl::init();

  //Decomp decomp;
  //decomp.initialize();

  std::cout << "Initializing mesh" << std::endl;
  Mesh mesh;
  mesh.read("initial_state.nc");

  State state("initial_state.nc", mesh);
  state.compute_ssh(mesh);
  
  MPI_Barrier(MPI_COMM_WORLD);
  io.create("output.nc", state);
  io.write(state);

  Timestep timestep(dt, mesh);
  
  yakl::fence();

  while (t < end_time) {

    std::cout << "Time: " << t << std::endl;
    timestep.RKStep(state, mesh, t);

    t = t + dt;
    state.compute_ssh(mesh);
    io.write(state);
  }

  yakl::finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}	
