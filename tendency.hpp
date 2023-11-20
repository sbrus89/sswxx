#include "const.h"
#include <stdlib.h>
#include <stdio.h>

class Tendency{
public:

  real2dHost layerThickness;
  real2dHost normalVelocity;
  real1dHost ssh;
  
  int nCells;
  int nEdges;
  int nVertLevels;
  
  Tendency(Mesh &mesh) {
  
    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;  
  
    layerThickness = real2dHost("layerThickness_tend", nCells, nVertLevels);
    normalVelocity = real2dHost("normalVelocity_tend", nEdges, nVertLevels);
    ssh = real1dHost("ssh", nCells);
  }
  
  void layerThicknessTendencies(Mesh &mesh, real2dHost &layerThickness_stage, real2dHost &normalVelocity_stage, real t) {
  
    parallel_for(SimpleBounds<2>(nCells, nVertLevels), YAKL_LAMBDA(int iCell, int kLevel) {
      
      real divergence = 0.0;
      for (int i=0; i<mesh.nEdgesOnCell(iCell); i++) {
        int iEdge = mesh.edgesOnCell(iCell,i);
        int cell2 = mesh.cellsOnCell(iCell,i);
        real hAvg = 0.5*(mesh.bottomDepth(iCell) + mesh.bottomDepth(cell2));
        divergence = divergence + mesh.edgeSignOnCell(iCell,i) * hAvg * normalVelocity_stage(iEdge,kLevel) * mesh.dvEdge(iEdge);
      }

      layerThickness(iCell,kLevel) = divergence/mesh.areaCell(iCell);
    });

  }
  
  void normalVelocityTendencies(Mesh &mesh, real2dHost &layerThickess_stage, real2dHost &normalVelocity_stage, real t){
  
    int iEdge, kLevel;
    int cell1, cell2;
    int i, j, eoe;

    parallel_for(SimpleBounds<1>(nCells), YAKL_LAMBDA(int iCell){

      real totalThickness = 0.0;
      for (int kLevel=0; kLevel<nVertLevels; kLevel++) {
        totalThickness = totalThickness + layerThickess_stage(iCell,kLevel);
      }
      ssh(iCell) = totalThickness - mesh.bottomDepth(iCell);
      
    });

    parallel_for(SimpleBounds<2>(nEdges, nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel){

      int cell1 = mesh.cellsOnEdge(iEdge,0); 
      int cell2 = mesh.cellsOnEdge(iEdge,1);

      normalVelocity(iEdge,kLevel) = -gravity*(ssh(cell2) - ssh(cell1))/mesh.dcEdge(iEdge);
    });
  
    parallel_for(SimpleBounds<2>(nEdges, nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel){

      real qe = 0.0;
      for (int j=0; j<mesh.nEdgesOnEdge(iEdge); j++) {
        int eoe = mesh.edgesOnEdge(iEdge,j);
        qe = qe + mesh.weightsOnEdge(iEdge,j)*normalVelocity_stage(eoe,kLevel)*mesh.fEdge(eoe);
      }

      normalVelocity(iEdge,kLevel) = normalVelocity(iEdge,kLevel) + qe;
    });

  }

};
