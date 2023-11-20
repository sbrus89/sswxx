#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "const.h"
#include "pnetcdf.h"
#include "mesh.hpp"
#include "state.hpp"
#include "timestep.hpp"
//#include "tendency.hpp"
#include <ctime>
#include <chrono>


int main (int argc, char **argv) {
  int nranks;
  int myrank;
  real dt, t;
  real end_time;

  end_time = 10.0*3600.0;
  dt = 20.0*60.0;
  t = 0.0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  std::cout << "Number of ranks: " << nranks << "\n";


  
  std::cout << "Initializing yakl\n";
  yakl::init();

  std::cout << "Initializing mesh \n";
  Mesh mesh;
  mesh.read("initial_state.nc");

  State state("initial_state.nc", mesh);
  
  MPI_Barrier(MPI_COMM_WORLD);

  Timestep timestep(dt, mesh);
  
  yakl::fence();

  while (t < end_time) {

    std::cout << "Time: " << t << "\n";
    timestep.RKStep(state, mesh, t);

    t = t + dt;
  }

  yakl::finalize();
  MPI_Finalize();
}	
