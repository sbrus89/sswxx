#include "const.h"
//#include "YAKL_netcdf.h"
#include <stdlib.h>
#include <stdio.h>

template<int M>
class State{
public:

  using int1d = yakl::Array<int   ,1, M>; 
  using int2d = yakl::Array<int   ,2, M>; 
  using int3d = yakl::Array<int   ,3, M>; 
  using real1d = yakl::Array<real  ,1, M>; 
  using real2d = yakl::Array<real  ,2, M>; 
  using real3d = yakl::Array<real  ,3, M>; 

  real2d layerThickness, layerThickness_new;
  real2d normalVelocity, normalVelocity_new;
  real1d ssh;

  int nCells;
  int nEdges;
  int nVertLevels;

  State (const char *mesh_file, Mesh<M> &mesh) {
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
    layerThickness = io.read<real2d>("layerThickness", __LINE__);
    normalVelocity = io.read<real2d>("normalVelocity", __LINE__);
    io.close();
    
    std::cout << "done with initial conditions" << std::endl;

    layerThickness_new = real2d("layerThickness_new", nCells, nVertLevels);
    normalVelocity_new = real2d("normalVelocity_new", nEdges, nVertLevels);
    ssh = real1d("ssh", nCells);
    

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

  void compute_ssh(Mesh<M> &mesh) {
  
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
