#include "const.h"
#include <stdlib.h>
#include <stdio.h>

template<int M>
class Mesh{
public:
   int nCells, nEdges, nVertices;
   int nVertLevels;
   int maxEdges, maxEdges2;
   int vertexDegree;

   using int1d = yakl::Array<int   ,1, M>; 
   using int2d = yakl::Array<int   ,2, M>; 
   using int3d = yakl::Array<int   ,3, M>; 
   using real1d = yakl::Array<real  ,1, M>; 
   using real2d = yakl::Array<real  ,2, M>; 
   using real3d = yakl::Array<real  ,3, M>; 

   int2d edgesOnEdge, cellsOnEdge, verticesOnEdge;
   int2d cellsOnVertex, edgesOnVertex;
   int2d cellsOnCell, edgesOnCell, verticesOnCell;
   int1d indexToEdgeID, indexToVertexID, indexToCellID;
   int1d nEdgesOnEdge, nEdgesOnCell;
   int1d boundaryVertex;
   int1d obtuseTriangle;
  
   real2d weightsOnEdge;
   real1d xCell, yCell, zCell;
   real1d lonCell, latCell;
   real1d xEdge, yEdge, zEdge;
   real1d lonEdge, latEdge;
   real1d xVertex, yVertex, zVertex;
   real1d lonVertex, latVertex;
   real1d dcEdge, dvEdge;
   real1d angleEdge;
   real1d cellQuality, triangleQuality, triangleAngleQuality;
   real1d areaCell, areaTriangle;
   real1d gridSpacing, meshDensity;
   real2d kiteAreasOnVertex;
   real1d bottomDepth;
   real1d fEdge, fVertex, fCell;

   real2d edgeSignOnCell;


   void read(const char* mesh_file){
     IO io;

     std::cout << mesh_file << std::endl;

     io.open(mesh_file);
     nCells = io.read_dim("nCells", __LINE__);
     nEdges = io.read_dim("nEdges", __LINE__);
     nVertices = io.read_dim("nVertices", __LINE__);
     nVertLevels= io.read_dim("nVertLevels", __LINE__);
     maxEdges = io.read_dim("maxEdges", __LINE__);
     maxEdges2 = io.read_dim("maxEdges2", __LINE__);
     vertexDegree = io.read_dim("vertexDegree", __LINE__);

     xCell = io.read<real1d>("xCell", __LINE__); 
     yCell = io.read<real1d>("yCell", __LINE__); 
     zCell = io.read<real1d>("zCell", __LINE__); 
     lonCell = io.read<real1d>("lonCell", __LINE__); 
     latCell = io.read<real1d>("latCell", __LINE__); 

     xEdge = io.read<real1d>("xEdge", __LINE__); 
     yEdge = io.read<real1d>("yEdge", __LINE__); 
     zEdge = io.read<real1d>("zEdge", __LINE__); 
     lonEdge = io.read<real1d>("lonEdge", __LINE__); 
     latEdge = io.read<real1d>("latEdge", __LINE__); 

     xVertex = io.read<real1d>("xVertex", __LINE__); 
     yVertex = io.read<real1d>("yVertex", __LINE__); 
     zVertex = io.read<real1d>("zVertex", __LINE__); 
     lonVertex = io.read<real1d>("lonVertex", __LINE__); 
     latVertex = io.read<real1d>("latVertex", __LINE__); 

     weightsOnEdge = io.read<real2d>("weightsOnEdge", __LINE__); 

     angleEdge = io.read<real1d>("angleEdge", __LINE__); 
     dcEdge = io.read<real1d>("dcEdge", __LINE__); 
     dvEdge = io.read<real1d>("dvEdge", __LINE__); 

     kiteAreasOnVertex = io.read<real2d>("kiteAreasOnVertex", __LINE__); 
     areaTriangle = io.read<real1d>("areaTriangle", __LINE__); 
     areaCell= io.read<real1d>("areaCell", __LINE__); 
     triangleQuality = io.read<real1d>("triangleQuality", __LINE__); 
     triangleAngleQuality = io.read<real1d>("triangleAngleQuality", __LINE__); 
     cellQuality = io.read<real1d>("cellQuality", __LINE__); 
     gridSpacing = io.read<real1d>("gridSpacing", __LINE__); 
     meshDensity = io.read<real1d>("meshDensity", __LINE__); 

     edgesOnEdge = io.read<int2d>("edgesOnEdge", __LINE__); 
     cellsOnEdge = io.read<int2d>("cellsOnEdge", __LINE__); 
     verticesOnEdge = io.read<int2d>("verticesOnEdge", __LINE__); 
     cellsOnVertex = io.read<int2d>("cellsOnVertex", __LINE__); 
     edgesOnVertex = io.read<int2d>("edgesOnVertex", __LINE__); 
     cellsOnCell = io.read<int2d>("cellsOnCell", __LINE__); 
     edgesOnCell = io.read<int2d>("edgesOnCell", __LINE__); 
     verticesOnCell = io.read<int2d>("verticesOnCell", __LINE__); 
     nEdgesOnEdge = io.read<int1d>("nEdgesOnEdge", __LINE__); 
     nEdgesOnCell = io.read<int1d>("nEdgesOnCell", __LINE__); 

     indexToEdgeID = io.read<int1d>("indexToEdgeID", __LINE__); 
     indexToCellID = io.read<int1d>("indexToCellID", __LINE__); 
     indexToVertexID = io.read<int1d>("indexToVertexID", __LINE__); 

     boundaryVertex = io.read<int1d>("boundaryVertex", __LINE__); 
     obtuseTriangle = io.read<int1d>("obtuseTriangle", __LINE__); 

     fEdge = io.read<real1d>("fEdge", __LINE__); 
     fCell= io.read<real1d>("fCell", __LINE__); 
     fVertex= io.read<real1d>("fVertex", __LINE__); 

     bottomDepth= io.read<real1d>("bottomDepth", __LINE__); 

     io.close();

     std::cout << "done reading mesh" << std::endl;

     fixIndices();
     computeEdgeSign();
     
  }

private:   

  void computeEdgeSign(void) {

    int iCell, jEdge;
    int nEdge;
    int j;

    edgeSignOnCell = real2d("edgeSignOnCell", nCells, maxEdges);
    
    for (iCell=0; iCell<nCells; iCell++) {
      nEdge = nEdgesOnCell(iCell);
      for (j=0; j<nEdge; j++) {
        jEdge = edgesOnCell(iCell,j);
        if (iCell == cellsOnEdge(jEdge,0)) {
          edgeSignOnCell(iCell,j) = -1.0;
        }
        else {
          edgeSignOnCell(iCell,j) = 1.0;
        }

      }
    }
  }

  void fixIndices(void) {

    int iCell, iEdge, iVertex;
    int j;

    for (iCell=0; iCell<nCells; iCell++) {
       for (j=0; j<maxEdges; j++) {
          cellsOnCell(iCell, j) = cellsOnCell(iCell, j) - 1;
          edgesOnCell(iCell, j) = edgesOnCell(iCell, j) - 1;
          verticesOnCell(iCell, j) = verticesOnCell(iCell, j) - 1;
       }
    }

   for (iEdge=0; iEdge<nEdges; iEdge++) {
      for(j=0; j<maxEdges2; j++) {
         edgesOnEdge(iEdge, j) = edgesOnEdge(iEdge, j) - 1;
      }
      for (j=0; j<2; j++) {
         cellsOnEdge(iEdge, j) = cellsOnEdge(iEdge, j) - 1;
         verticesOnEdge(iEdge, j) = verticesOnEdge(iEdge, j) - 1;
      }
   }

   for (iVertex=0; iVertex<nVertices; iVertex++) {
      for (j=0; j<vertexDegree; j++) {
         cellsOnVertex(iVertex, j) = cellsOnVertex(iVertex, j) - 1;
         edgesOnVertex(iVertex, j) = edgesOnVertex(iVertex, j) - 1;
      }
   }
  }


};
