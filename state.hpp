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

  State (Mesh<M> &mesh) {

    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;

    layerThickness = real2d("layerThickness", nCells, nVertLevels);
    normalVelocity = real2d("normalVelocity", nEdges, nVertLevels);
    layerThickness_new = real2d("layerThickness_new", nCells, nVertLevels);
    normalVelocity_new = real2d("normalVelocity_new", nEdges, nVertLevels);
    ssh = real1d("ssh", nCells);

  }

  void read_initial_condition(const char *mesh_file) {
    IO io;

    std::cout << "Reading initial conditions" << std::endl;
    std::cout << "  nCells: " << nCells << std::endl;
    std::cout << "  nEdges: " << nEdges << std::endl;
    std::cout << "  nVertLevels: " << nVertLevels << std::endl;

    std::cout << "begin reading initial conditions" << std::endl;

    io.open(mesh_file);
/*
    layerThickness = io.read<real2d>("layerThickness", __LINE__);
    normalVelocity = io.read<real2d>("normalVelocity", __LINE__);
*/
    io.read<real2d>("layerThickness", layerThickness, __LINE__);
    io.read<real2d>("normalVelocity", normalVelocity, __LINE__);
    io.close();
    
    std::cout << "done with initial conditions" << std::endl;

  }

  void update_time_level() {

    YAKL_SCOPE(layerThickness, this->layerThickness);
    YAKL_SCOPE(layerThickness_new, this->layerThickness_new);
    YAKL_SCOPE(nCells, this->nCells);
    YAKL_SCOPE(nVertLevels, this->nVertLevels);
    parallel_for(SimpleBounds<2>(nCells, nVertLevels), YAKL_LAMBDA(int iCell, int kLevel) {
        layerThickness(iCell,kLevel) = layerThickness_new(iCell,kLevel);
        layerThickness_new(iCell,kLevel) = 0.0;
    });

    YAKL_SCOPE(normalVelocity, this->normalVelocity);
    YAKL_SCOPE(normalVelocity_new, this->normalVelocity_new);
    YAKL_SCOPE(nEdges, this->nEdges);
    parallel_for(SimpleBounds<2>(nEdges, nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel) {
        normalVelocity(iEdge,kLevel) = normalVelocity_new(iEdge,kLevel);
        normalVelocity_new(iEdge,kLevel) = 0.0;
    });

  }

  void compute_ssh(Mesh<M> &mesh) {
  
    YAKL_SCOPE(ssh, this->ssh);
    YAKL_SCOPE(layerThickness, this->layerThickness);
    YAKL_SCOPE(nCells, this->nCells);
    YAKL_SCOPE(nVertLevels, this->nVertLevels);
    parallel_for(SimpleBounds<1>(nCells), YAKL_LAMBDA(int iCell) {
      real totalThickness = 0.0;
      for (int kLevel=0; kLevel<nVertLevels; kLevel++) {
        totalThickness = totalThickness + layerThickness(iCell,kLevel);
      }
      ssh(iCell) = totalThickness - mesh.bottomDepth(iCell);
    });

  }

  template<int N>
  void copy_to(State<N> state){

    state.nCells = nCells;
    state.nEdges = nEdges;
    state.nVertLevels = nVertLevels;

    layerThickness.deep_copy_to(state.layerThickness);
    layerThickness_new.deep_copy_to(state.layerThickness_new);
    normalVelocity.deep_copy_to(state.normalVelocity);
    normalVelocity_new.deep_copy_to(state.normalVelocity_new);
    ssh.deep_copy_to(state.ssh);
    
  }

};
