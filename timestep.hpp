#include "const.h"
#include <stdlib.h>
#include <stdio.h>
//#include "state.hpp"
//#include "mesh.hpp"
#include "tendency.hpp"

class Timestep{
public:

  real1dHost a, b, c;
  double dt;
  Tendency tend;
  int nStage;
  real2dHost layerThickness_stage;
  real2dHost normalVelocity_stage;

  Timestep (double dt, Mesh &mesh): tend(mesh) {

    nStage = 4;
     
    a = real1dHost("a",nStage);
    b = real1dHost("b",nStage);
    c = real1dHost("c",nStage);

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

    this->dt = dt;

    layerThickness_stage = real2dHost("layerThickness_stage", mesh.nCells, mesh.nVertLevels);
    normalVelocity_stage = real2dHost("normalVelocity_stage", mesh.nEdges, mesh.nVertLevels);

  }

  void RKStep (State &state, Mesh &mesh, real t) {

    int iStage;

    for (iStage=0; iStage<nStage; iStage++) {
      compute_stage(iStage, t, state, mesh);
    }

    state.update_time_level();

  }

  void compute_stage(int stage, real t, State &state, Mesh &mesh) {

    int iCell;
    int iEdge;
    int jLevel;

    for (iCell=0; iCell<mesh.nCells; iCell++) {
      for (jLevel=0; jLevel<mesh.nVertLevels; jLevel++) {
        layerThickness_stage(iCell,jLevel) = state.layerThickness(iCell,jLevel) + a(stage)*dt*tend.layerThickness(iCell,jLevel);
      }
    }
    for (iEdge=0; iEdge<mesh.nEdges; iEdge++) {
      for (jLevel=0; jLevel<mesh.nVertLevels; jLevel++) {
        normalVelocity_stage(iEdge,jLevel) = state.normalVelocity(iEdge,jLevel) + a(stage)*dt*tend.normalVelocity(iEdge,jLevel);
      }
    }

    compute_tendencies(mesh, t+c(stage)*dt);
    accumulate_stage(stage, state);

  }

  void compute_tendencies(Mesh &mesh, double t) {

   tend.layerThicknessTendencies(mesh, layerThickness_stage, normalVelocity_stage, t);
   tend.normalVelocityTendencies(mesh, layerThickness_stage, normalVelocity_stage, t);

  }

  void accumulate_stage(int stage, State &state) {

    int iCell;
    int iEdge;
    int jLevel;

    if (stage==0) {
      for (iCell=0; iCell<state.nCells; iCell++) {
        for (jLevel=0; jLevel<state.nVertLevels; jLevel++) { 
          state.layerThickness_new(iCell,jLevel) = state.layerThickness(iCell,jLevel) + b(stage)*dt*tend.layerThickness(iCell,jLevel);
        }
      }
      for (iEdge=0; iEdge<state.nEdges; iEdge++) {
        for (jLevel=0; jLevel<state.nVertLevels; jLevel++) { 
          state.normalVelocity_new(iEdge,jLevel) = state.normalVelocity(iEdge,jLevel) + b(stage)*dt*tend.normalVelocity(iEdge,jLevel);
        }
      }

    }
    else {
      for (iCell=0; iCell<state.nCells; iCell++) {
        for (jLevel=0; jLevel<state.nVertLevels; jLevel++) { 
          state.layerThickness_new(iCell,jLevel) = state.layerThickness_new(iCell,jLevel) + b(stage)*dt*tend.layerThickness(iCell,jLevel);
        }
      }
      for (iEdge=0; iEdge<state.nEdges; iEdge++) {
        for (jLevel=0; jLevel<state.nVertLevels; jLevel++) { 
          state.normalVelocity_new(iEdge,jLevel) = state.normalVelocity_new(iEdge,jLevel) + b(stage)*dt*tend.normalVelocity(iEdge,jLevel);
        }
      }
  }
  }

};	
