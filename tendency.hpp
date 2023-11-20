#include "const.h"
#include <stdlib.h>
#include <stdio.h>

class Tendency{
public:

  real2dHost layerThickness;
  real2dHost normalVelocity;
  
  int nCells;
  int nEdges;
  int nVertLevels;
  
  Tendency(Mesh &mesh) {
  
    this->nCells = mesh.nCells;
    this->nEdges = mesh.nEdges;
    this->nVertLevels = mesh.nVertLevels;  
  
    layerThickness = real2dHost("layerThickness_tend", nCells, nVertLevels);
    normalVelocity = real2dHost("normalVelocity_tend", nEdges, nVertLevels);
  }
  
  void layerThicknessTendencies(Mesh &mesh, real2dHost &layerThickness_stage, real2dHost &normalVelocity_stage, real t) {
  
    int iCell, kLevel;
    int i, iEdge;
    int cell2;
    real hAvg;

    for (iCell=0; iCell<nCells; iCell++) {
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        layerThickness(iCell,kLevel) = 0.0;
      }
    }
  
    for (iCell=0; iCell<nCells; iCell++) {
      for (i=0; i<mesh.nEdgesOnCell(iCell); i++) {
        iEdge = mesh.edgesOnCell(iCell,i);
        cell2 = mesh.cellsOnCell(iCell,i);
        //hAvg = 0.5*(layerThickness_stage(iCell,kLevel) + layerThickness_stage(cell2,kLevel));
        hAvg = 0.5*(mesh.bottomDepth(iCell) + mesh.bottomDepth(cell2));
        for (kLevel=0; kLevel<nVertLevels; kLevel++) {
          layerThickness(iCell,kLevel) = layerThickness(iCell,kLevel) + mesh.edgeSignOnCell(iCell,i)*hAvg*normalVelocity_stage(iEdge,kLevel)*mesh.dvEdge(iEdge)/mesh.areaCell(iCell);
        }
      }
    }

  }
  
  void normalVelocityTendencies(Mesh &mesh, real2dHost &layerThickess_stage, real2dHost &normalVelocity_stage, real t){
  
    int iEdge, kLevel;
    int cell1, cell2;
    int i, j, eoe;
    real totalThickness1, totalThickness2;
    real ssh1, ssh2;

    for (iEdge=0; iEdge<nEdges; iEdge++) {
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        normalVelocity(iEdge,kLevel) = 0.0;
      }
    }
    
    for (iEdge=0; iEdge<nEdges; iEdge++) {

      cell1 = mesh.cellsOnEdge(iEdge,0); 
      cell2 = mesh.cellsOnEdge(iEdge,1);

      totalThickness1 = 0.0;
      totalThickness2 = 0.0;
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        totalThickness1 = totalThickness1 + layerThickess_stage(cell1,kLevel);
        totalThickness2 = totalThickness2 + layerThickess_stage(cell2,kLevel);
      }
      ssh1 = totalThickness1 - mesh.bottomDepth(cell1);
      ssh2 = totalThickness2 - mesh.bottomDepth(cell2);
      
      
      for (kLevel=0; kLevel<nVertLevels; kLevel++) {
        normalVelocity(iEdge,kLevel) = normalVelocity(iEdge,kLevel) - gravity*(ssh2-ssh1)/mesh.dcEdge(iEdge);
      }
    }
  
    for (iEdge=0; iEdge<nEdges; iEdge++) {
      for (j=0; j<mesh.nEdgesOnEdge(iEdge); j++) {
        eoe = mesh.edgesOnEdge(iEdge,j);
        for (kLevel=0; kLevel<nVertLevels; kLevel++) {
          normalVelocity(iEdge,kLevel) = normalVelocity(iEdge,kLevel) + mesh.weightsOnEdge(iEdge,j)*normalVelocity_stage(eoe,kLevel)*mesh.fEdge(eoe);
        }
      }
    }

  }

};
