#include "const.h"
#include <stdlib.h>
#include <stdio.h>
//#include "state.hpp"
//#include "mesh.hpp"
#include "tendency.hpp"

template<int M>
class Timestep{
public:

  using int1d = yakl::Array<int   ,1, M>; 
  using int2d = yakl::Array<int   ,2, M>; 
  using int3d = yakl::Array<int   ,3, M>; 
  using real1d = yakl::Array<real  ,1, M>; 
  using real2d = yakl::Array<real  ,2, M>; 
  using real3d = yakl::Array<real  ,3, M>; 

  real a[4], b[4], c[4];
  real dt;
  Tendency<M> tend;
  int nStage;
  real2d layerThickness_stage;
  real2d normalVelocity_stage;

  Timestep (real dt, Mesh<M> &mesh): tend(mesh) {

    nStage = 4;

    a[0] = 0.0;
    a[1] = 1.0/2.0;
    a[2] = 1.0/2.0;
    a[3] = 1.0;

    b[0] = 1.0/6.0;
    b[1] = 1.0/3.0;
    b[2] = 1.0/3.0;
    b[3] = 1.0/6.0;

    c[0] = 0.0;
    c[1] = 1.0/2.0;
    c[2] = 1.0/2.0;
    c[3] = 1.0;

    this->dt = dt;

    layerThickness_stage = real2d("layerThickness_stage", mesh.nCells, mesh.nVertLevels);
    normalVelocity_stage = real2d("normalVelocity_stage", mesh.nEdges, mesh.nVertLevels);

  }

  void RKStep (State<M> &state, Mesh<M> &mesh, real t) {

    int iStage;

    for (iStage=0; iStage<nStage; iStage++) {
      compute_stage(iStage, t, a[iStage], b[iStage], c[iStage], dt, state, mesh);
    }

    state.update_time_level();

  }

  void compute_stage(int stage, real t, real a_stage, real b_stage, real c_stage, real dt,  State<M> &state, Mesh<M> &mesh) {
    
    YAKL_SCOPE(layerThickness_stage, this->layerThickness_stage);
    YAKL_SCOPE(layerThickness_tend, tend.layerThickness);
    parallel_for(SimpleBounds<2>(mesh.nCells,mesh.nVertLevels), YAKL_LAMBDA(int iCell, int kLevel){
      layerThickness_stage(iCell,kLevel) = state.layerThickness(iCell,kLevel)+ a_stage*dt*layerThickness_tend(iCell,kLevel);
    });

    YAKL_SCOPE(normalVelocity_stage, this->normalVelocity_stage);
    YAKL_SCOPE(normalVelocity_tend, tend.normalVelocity);
    parallel_for(SimpleBounds<2>(mesh.nEdges,mesh.nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel){
        normalVelocity_stage(iEdge,kLevel) = state.normalVelocity(iEdge,kLevel) + a_stage*dt*normalVelocity_tend(iEdge,kLevel);
    });

    compute_tendencies(mesh, t+c_stage*dt);
    accumulate_stage(stage, b_stage, dt, state);

  }

  void compute_tendencies(Mesh<M> &mesh, double t) {

   tend.layerThicknessTendencies(mesh, layerThickness_stage, normalVelocity_stage, t);
   tend.normalVelocityTendencies(mesh, layerThickness_stage, normalVelocity_stage, t);

  }

  void accumulate_stage(int stage, real b_stage, real dt, State<M> &state) {

    if (stage==0) {

      YAKL_SCOPE(layerThickness_tend, tend.layerThickness);
      parallel_for(SimpleBounds<2>(state.nCells,state.nVertLevels), YAKL_LAMBDA(int iCell, int kLevel){
        state.layerThickness_new(iCell,kLevel) = state.layerThickness(iCell,kLevel) + b_stage*dt*layerThickness_tend(iCell,kLevel);
      });

      YAKL_SCOPE(normalVelocity_tend, tend.normalVelocity);
      parallel_for(SimpleBounds<2>(state.nEdges,state.nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel){
        state.normalVelocity_new(iEdge,kLevel) = state.normalVelocity(iEdge,kLevel) + b_stage*dt*normalVelocity_tend(iEdge,kLevel);
      });

    }
    else {

      YAKL_SCOPE(layerThickness_tend, tend.layerThickness);
      parallel_for(SimpleBounds<2>(state.nCells,state.nVertLevels), YAKL_LAMBDA(int iCell, int kLevel){
        state.layerThickness_new(iCell,kLevel) = state.layerThickness_new(iCell,kLevel) + b_stage*dt*layerThickness_tend(iCell,kLevel);
      });

      YAKL_SCOPE(normalVelocity_tend, tend.normalVelocity);
      parallel_for(SimpleBounds<2>(state.nEdges,state.nVertLevels), YAKL_LAMBDA(int iEdge, int kLevel){
        state.normalVelocity_new(iEdge,kLevel) = state.normalVelocity_new(iEdge,kLevel) + b_stage*dt*normalVelocity_tend(iEdge,kLevel);
      });

    }
  }

};	
