#include "const.h"
#include "YAKL_netcdf.h"
#include <stdlib.h>
#include <stdio.h>

class State{
public:

  real2dHost layerThickness, layerThickness_new;
  real2dHost normalVelocity, normalVelocity_new;

  int nCells;
  int nEdges;
  int nVertLevels;

  State (const char *mesh_file, Mesh &mesh) {
    IO io;
    //yakl::SimpleNetCDF nc;

    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;

    std::cout << "Reading initial conditions\n";
    std::cout << "  nCells: " << nCells << "\n";
    std::cout << "  nEdges: " << nEdges << "\n";
    std::cout << "  nVertLevels: " << nVertLevels << "\n";


    std::cout << "begin reading initial conditions\n";

    io.open(mesh_file);
    layerThickness = io.read<real2dHost>("layerThickness", __LINE__);
    normalVelocity = io.read<real2dHost>("normalVelocity", __LINE__);
    io.close();
    
    //nc.open(mesh_file);
    //nc.read(layerThickness, "layerThickness");
    //nc.read(normalVelocity, "normalVelocity");
    //nc.close();


    std::cout << "done with initial conditions\n";

    layerThickness_new = real2dHost("layerThickness_new", nCells, nVertLevels);
    normalVelocity_new = real2dHost("normalVelocity_new", nEdges, nVertLevels);
    

  }

  void halo_exchange() {

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

};	
