#include "const.h"
//#include "YAKL_netcdf.h"
#include <stdlib.h>
#include <stdio.h>

class State{
public:

  real2dHost layerThickness, layerThickness_new;
  real2dHost normalVelocity, normalVelocity_new;
  real1dHost ssh;

  int nCells;
  int nEdges;
  int nVertLevels;

  State (const char *mesh_file, Mesh &mesh) {
    IO io;

    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;

    std::cout << "Reading initial conditions" << std::endl;
    std::cout << "  nCells: " << nCells << std::endl;
    std::cout << "  nEdges: " << nEdges << std::endl;
    std::cout << "  nVertLevels: " << nVertLevels << std::endl;


    std::cout << "begin reading initial conditions" << std::endl;

    io.open(mesh_file);
    layerThickness = io.read<real2dHost>("layerThickness", __LINE__);
    normalVelocity = io.read<real2dHost>("normalVelocity", __LINE__);
    io.close();
    
    std::cout << "done with initial conditions" << std::endl;

    layerThickness_new = real2dHost("layerThickness_new", nCells, nVertLevels);
    normalVelocity_new = real2dHost("normalVelocity_new", nEdges, nVertLevels);
    ssh = real1dHost("ssh", nCells);
    

  }

  void update_time_level() {

    int iCell;
    int iEdge;
    int kLevel;

    for (iCell=0; iCell<nCells; iCell++) {
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        layerThickness(iCell,kLevel) = layerThickness_new(iCell,kLevel);
        layerThickness_new(iCell,kLevel) = 0.0;
      }
    }

    for (iEdge=0; iEdge<nEdges; iEdge++) {
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        normalVelocity(iEdge,kLevel) = normalVelocity_new(iEdge,kLevel);
        normalVelocity_new(iEdge,kLevel) = 0.0;
      }
    }

  }

  void compute_ssh(Mesh &mesh) {
  
  int iCell;
  int kLevel;
  real totalThickness;

  for (iCell=0; iCell<nCells; iCell++) {
    totalThickness = 0.0;
    for (kLevel=0; kLevel<nVertLevels; kLevel++) {
      totalThickness = totalThickness + layerThickness(iCell,kLevel);
    }
    ssh(iCell) = totalThickness - mesh.bottomDepth(iCell);
  }

  }

};	
