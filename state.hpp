#include "const.h"
#include "io.hpp"
#include <stdlib.h>
#include <stdio.h>

class State{
public:

  real3d ssh;
  real3d normalVelocity;
  real3d temp;
  real3d dens;
  IO io;

  int nCells;
  int nEdges;
  int nVertLevels;

  void State(int, nTimeLevels, int nCells, int nEdges, int nVertLevels) {

    ssh = real3d("ssh", nTimeLevels, nCells, nVertLevels);
    normalVelocity = real3D("normalVelocity", nTimeLevels, nEdges, nVertLevels);
    temp = real3d("temp", nTimeLevels, nCells, nVertLevels);
    dens = real3d("dens", nTimeLevels, nCells, nVertLevels);

    this->nCells = nCells;
    this->nEdges = nEdges;
    this->nVertLevels = nVertLevels;

  }

  void halo_exchange() {

  }

  void update_time_level() {

    int iCell;
    int iEdge;
    int jLevel;

    for (iCell=0; iCell<nCells; iCell++) {
      for (jLevel=0; jLevel<nVertLevels; jLevel++) {
        ssh(1,iCell,jLevel) = ssh(2,iCell,jLevel)
      }
    }
    for (iEdge=0; iEdge<nEdge; iEdge++) {
      for (jLevel=0; jLevels<nVertLevels; jLevel++) {
        normalVelocity(1,iCell,jLevel) = normalVelocity(2,iCell,jLevel)
      }
    }

  }

  real2dHost read_double_var(const char* var_str, int dim1, int dim2, int dim3, int line) {
    int var_id;
    int i, j, k, l;
    real2dHost var;

    ncwrap(ncmpi_inq_varid(nc_id, var_str, &var_id), line);
    double buff[dim1*dim2];
    ncwrap(ncmpi_get_var_double_all(nc_id, var_id, buff), line);
    var = real2dHost(var_str, dim1, dim2);

    std::cout << "\n";
    std::cout << var_str << "\n";

    l = 0;
    for (j=0; j<dim1; j++) {
      for (i=0; i<dim2; i++) {
        for (k=0; k<dim3; k++) {
          var(j,i,k) = buff[l];
          std::cout << var(j,i,k) << ", ";
          l++;
        }
      }   
      std::cout << "\n";
    }   
    return var;
  }

};	
