#include "const.h"
#include <stdlib.h>
#include <stdio.h>
#include "state.hpp"
#include "mesh.hpp"

class Timestep{
public:

  real1d a, b, c
  State state, tend; 
  double dt;
  int tNew, tOld;

  void Timestep (double dt, Mesh &mesh): tend(1, mesh.nCells, mesh.nEdges, mesh.nVertLevels),
                                         state(2, mesh.nCells, mesh.nEdges, mesh.nVertLevels)

{

    nStage = 4;
     
    a = real1d('a',nStage);
    b = real1d('b',nStage);
    c = real1d('c',nStage);

    a(0) = 0.0;
    a(1) = 1.0/2.0;
    a(2) = 1.0/2.0;
    a(3) = 1.0;

    b(0) = 1.0/6.0;
    b(1) = 1.0/3.0;
    b(2) = 1.0/3.0;
    b(3) = 1.0/6.0;

    c(0) = 0.0;
    c(1) = 1.0/2.0;
    c(2) = 1.0/2.0;
    c(3) = 1.0;

    //tend = new State(2, mesh.nCells, mesh.nEdges, mesh.nVertLevels);
    this->dt = dt;

    tNew = 1;
    tOld = 0;

    }

  void RKStep (State &state, Mesh &mesh, t) {

    int iStage;

    for (iStage=0; i<nStage; iStage++) {
      compute_substep(iStage, state, mesh, tend);
    }

  }

  void compute_substep(int stage, State &state, Mesh &mesh, State &tend) {

    int iCell;
    int jLevel;

    for (iCell=0; iCell<mesh.nCells; iCell++) {
      for (jLevel=0; jLevel<mesh.nVertLevels, jLevel++) {
        tend.ssh(1,iCell,jLevel) = state.ssh(1,iCell,jLevel) + a(stage)*tend.ssh(2,iCell,jLevel);
        tend.ssh(2,iCell,jLevel) = 0.0;
      }
    }
   for (iEdge=0; iEdge<mesh.nEdges; iEdge++) {
     for (jLevel=0; jLevel<mesh.nVertLevels, jLevel++) {
       tend.normalVelocity(1,iEdge,jLevel) = state.normalVelocity(1,iCell,jLevel) + a(stage)*tend.normalVelocity(2,iEdge,jLevel);
       tend.normalVelocity(2,iEdge,jLevel) = 0.0;
     }
   }

    compute_tendencies(tend, mesh, t+c(stage)*dt);
    accumulate_substep(stage, state, tend);
  }

  void compute_tendencies(State &tend, Mesh &mesh, double t) {

   int iCell;
   int j, jEdge;
   int kLevel;

   for (iCell=0; iCell<mesh.nCells; iCell++) {
     nEdges = mesh.nEdgesOnCell(1,iCell);
     for (j=0; j<nEdges; j++) {
       jEdge = mesh.edgesOnCell(j,iCell);
       for (kLevel=0; kLevel<mesh.nVertLevels; kLevel++) {
         tend.ssh(2,iCell,jLevel) = tend.ssh(2,iCell,jLevel) + 
       }
     }
   } 

  }

  void accumulate_subtetps(int stage, State &state, State &tend) {

    int iCell;
    int iEdge;
    int jLevel;

    for (iCell=0; i<state.nCells; iCell++) {
      for (jLevel=0; j<state.nVertLevels; jLevel++) { 
        state.ssh(2,iCell,jLevel) = state.ssh(2,iCell,jLevel) + b(stage)*dt*tend.ssh(2,iCell,jLevel);
      }
    }
    for (iEdge=0; i<state.nEdges; iEdge++) {
      for (jLevel=0; j<state.nVertLevels; jLevel++) { 
        state.normalVelocity(2,iEdge,jLevel) = state.normalVelocity(2,iEdge,jLevel) + b(stage)*dt*tend.normalVelocity(2,iEdge,jLevel);
      }
    }

  }

};	
