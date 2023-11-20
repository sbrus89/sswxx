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
//#include "decomp.hpp"
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
  Mesh<yakl::memHost> mesh_host;
  mesh_host.read("initial_state.nc");

  Mesh<yakl::memDevice> mesh;
  mesh_host.copy_to_device(mesh);

  State<yakl::memHost> state_host(mesh_host);
  state_host.read_initial_condition("initial_state.nc");

  std::cout << "Initializing device state" << std::endl;
  State<yakl::memDevice> state(mesh);
  state_host.copy_to(state);
  state.compute_ssh(mesh);
  state.copy_to(state_host);

  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "Creating output file" << std::endl;
  io.create("output.nc", state_host);
  std::cout << "writing initial condition" << std::endl;
  io.write(state_host);

  std::cout << "initializing timestep" << std::endl;
  Timestep<yakl::memDevice> timestep(dt, mesh);
  std::cout << "done" << std::endl;
  
  yakl::fence();

  while (t < end_time) {

    std::cout << "Time: " << t << std::endl;
    timestep.RKStep(state, mesh, t);
    yakl::fence();

    t = t + dt;
    state.compute_ssh(mesh);
    state.copy_to(state_host);
    io.write(state_host);
    yakl::fence();
  }

  yakl::finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}	
