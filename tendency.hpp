#include "const.h"
#include <stdlib.h>
#include <stdio.h>

template<int M>
class Tendency{
public:

  using int1d = yakl::Array<int   ,1, M>; 
  using int2d = yakl::Array<int   ,2, M>; 
  using int3d = yakl::Array<int   ,3, M>; 
  using real1d = yakl::Array<real  ,1, M>; 
  using real2d = yakl::Array<real  ,2, M>; 
  using real3d = yakl::Array<real  ,3, M>; 

  real2d layerThickness;
  real2d normalVelocity;
  real1d ssh;
  
  int nCells;
  int nEdges;
  int nVertLevels;
  
  Tendency(Mesh<M> &mesh) {
  
    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;  
  
    layerThickness = real2d("layerThickness_tend", nCells, nVertLevels);
    normalVelocity = real2d("normalVelocity_tend", nEdges, nVertLevels);
    ssh = real1d("ssh", nCells);
  }
  
  void layerThicknessTendencies(Mesh<M> &mesh, real2d &layerThickness_stage, real2d &normalVelocity_stage, real t) {
  
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
  
  void normalVelocityTendencies(Mesh<M> &mesh, real2d &layerThickess_stage, real2d &normalVelocity_stage, real t){
  
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
