#include "const.h"
#include "io.hpp"
#include <stdlib.h>
#include <stdio.h>

class Mesh{
public:
   int nCells, nEdges, nVertices;
   int maxEdges, maxEdges2;
   int vertexDegree;

   int2dHost edgesOnEdge, cellsOnEdge, verticesOnEdge;
   int2dHost cellsOnVertex, edgesOnVertex;
   int2dHost cellsOnCell, edgesOnCell, verticesOnCell;
   int1dHost indexToEdgeID, indexToVertexID, indexToCellID;
   int1dHost nEdgesOnEdge, nEdgesOnCell;
   int1dHost boundaryVertex;
   int1dHost obtuseTriangle;
  
   real2dHost weightsOnEdge;
   real1dHost xCell, yCell, zCell;
   real1dHost lonCell, latCell;
   real1dHost xEdge, yEdge, zEdge;
   real1dHost lonEdge, latEdge;
   real1dHost xVertex, yVertex, zVertex;
   real1dHost lonVertex, latVertex;
   real1dHost dcEdge, dvEdge;
   real1dHost angleEdge;
   real1dHost cellQuality, triangleQuality, triangleAngleQuality;
   real1dHost areaCell, areaTriangle;
   real1dHost gridSpacing, meshDensity;
   real2dHost kiteAreasOnVertex;

   real2dHost edgeSignOnCell;


   void read(const char* mesh_file){
     IO io;
     std::cout << mesh_file << "\n";

     io.open(mesh_file);
     nCells = io.read_dim("nCells", __LINE__);
     nEdges = io.read_dim("nEdges", __LINE__);
     nVertices = io.read_dim("nVertices", __LINE__);
     maxEdges = io.read_dim("maxEdges", __LINE__);
     maxEdges2 = io.read_dim("maxEdges2", __LINE__);
     vertexDegree = io.read_dim("vertexDegree", __LINE__);

     xCell = io.read<real1dHost>("xCell", __LINE__); 
     yCell = io.read<real1dHost>("yCell", __LINE__); 
     zCell = io.read<real1dHost>("zCell", __LINE__); 
     lonCell = io.read<real1dHost>("lonCell", __LINE__); 
     latCell = io.read<real1dHost>("latCell", __LINE__); 

     xEdge = io.read<real1dHost>("xEdge", __LINE__); 
     yEdge = io.read<real1dHost>("yEdge", __LINE__); 
     zEdge = io.read<real1dHost>("zEdge", __LINE__); 
     lonEdge = io.read<real1dHost>("lonEdge", __LINE__); 
     latEdge = io.read<real1dHost>("latEdge", __LINE__); 

     xVertex = io.read<real1dHost>("xVertex", __LINE__); 
     yVertex = io.read<real1dHost>("yVertex", __LINE__); 
     zVertex = io.read<real1dHost>("zVertex", __LINE__); 
     lonVertex = io.read<real1dHost>("lonVertex", __LINE__); 
     latVertex = io.read<real1dHost>("latVertex", __LINE__); 

     weightsOnEdge = io.read<real2dHost>("weightsOnEdge", __LINE__); 

     angleEdge = io.read<real1dHost>("angleEdge", __LINE__); 
     dcEdge = io.read<real1dHost>("dcEdge", __LINE__); 
     dvEdge = io.read<real1dHost>("dvEdge", __LINE__); 

     kiteAreasOnVertex = io.read<real2dHost>("kiteAreasOnVertex", __LINE__); 
     areaTriangle = io.read<real1dHost>("areaTriangle", __LINE__); 
     areaCell= io.read<real1dHost>("areaCell", __LINE__); 
     triangleQuality = io.read<real1dHost>("triangleQuality", __LINE__); 
     triangleAngleQuality = io.read<real1dHost>("triangleAngleQuality", __LINE__); 
     cellQuality = io.read<real1dHost>("cellQuality", __LINE__); 
     gridSpacing = io.read<real1dHost>("gridSpacing", __LINE__); 
     meshDensity = io.read<real1dHost>("meshDensity", __LINE__); 

     edgesOnEdge = io.read<int2dHost>("edgesOnEdge", __LINE__); 
     cellsOnEdge = io.read<int2dHost>("cellsOnEdge", __LINE__); 
     verticesOnEdge = io.read<int2dHost>("verticesOnEdge", __LINE__); 
     cellsOnVertex = io.read<int2dHost>("cellsOnVertex", __LINE__); 
     edgesOnVertex = io.read<int2dHost>("edgesOnVertex", __LINE__); 
     cellsOnCell = io.read<int2dHost>("cellsOnCell", __LINE__); 
     edgesOnCell = io.read<int2dHost>("edgesOnCell", __LINE__); 
     verticesOnCell = io.read<int2dHost>("verticesOnCell", __LINE__); 
     nEdgesOnEdge = io.read<int1dHost>("nEdgesOnEdge", __LINE__); 
     nEdgesOnCell = io.read<int1dHost>("nEdgesOnCell", __LINE__); 

     indexToEdgeID = io.read<int1dHost>("indexToEdgeID", __LINE__); 
     indexToCellID = io.read<int1dHost>("indexToCellID", __LINE__); 
     indexToVertexID = io.read<int1dHost>("indexToVertexID", __LINE__); 

     boundaryVertex = io.read<int1dHost>("boundaryVertex", __LINE__); 
     obtuseTriangle = io.read<int1dHost>("obtuseTriangle", __LINE__); 

     edgeSignOnCell = computeEdgeSign();
     
     std::cout << "done reading mesh\n";
  }

private:   

  real2dHost computeEdgeSign(void) {

    int iCell, jEdge;
    int nEdge;
    int j;

    edgeSignOnCell = real2dHost("edgeSignOnCell", maxEdges, nCells);
    
    for (iCell=0; iCell<nCells; iCell++) {
      nEdge = nEdgesOnCell(iCell);
      for (j=0; j<nEdge; j++) {
        jEdge = edgesOnCell(j,iCell);
        if (iCell == cellsOnEdge(1,jEdge)) {
          edgeSignOnCell(j,iCell) = -1.0;
        }
        else {
          edgeSignOnCell(j,iCell) = 1.0;
        }

      }
    }

    return edgeSignOnCell;

  }

};	
