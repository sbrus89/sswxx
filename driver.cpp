#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "const.h"
#include "pnetcdf.h"
#include "mesh.hpp"
#include <ctime>
#include <chrono>


int main (int argc, char **argv) {
  int nranks;
  int myrank;
  int i;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  std::cout << "Number of ranks: " << nranks << "\n";

  yakl::init();

  Mesh mesh;
  mesh.read("mesh.nc");

  yakl::finalize();
  MPI_Finalize();
}	
