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

     std::cout << "begin reading arrays" << std::endl;
     io.read<real1d>("xCell", xCell, __LINE__); 
     io.read<real1d>("yCell", yCell, __LINE__); 
     io.read<real1d>("zCell", zCell, __LINE__); 
     io.read<real1d>("lonCell", lonCell, __LINE__); 
     io.read<real1d>("latCell", latCell, __LINE__); 

     io.read<real1d>("xEdge", xEdge, __LINE__); 
     io.read<real1d>("yEdge", yEdge, __LINE__); 
     io.read<real1d>("zEdge", zEdge, __LINE__); 
     io.read<real1d>("lonEdge", lonEdge, __LINE__); 
     io.read<real1d>("latEdge", latEdge, __LINE__); 

     io.read<real1d>("xVertex", xVertex, __LINE__); 
     io.read<real1d>("yVertex", yVertex, __LINE__); 
     io.read<real1d>("zVertex", zVertex, __LINE__); 
     io.read<real1d>("lonVertex", lonVertex, __LINE__); 
     io.read<real1d>("latVertex", latVertex, __LINE__); 

     io.read<real2d>("weightsOnEdge", weightsOnEdge, __LINE__); 

     io.read<real1d>("angleEdge", angleEdge, __LINE__); 
     io.read<real1d>("dcEdge", dcEdge, __LINE__); 
     io.read<real1d>("dvEdge", dvEdge, __LINE__); 

     io.read<real2d>("kiteAreasOnVertex", kiteAreasOnVertex, __LINE__); 
     io.read<real1d>("areaTriangle", areaTriangle, __LINE__); 
     io.read<real1d>("areaCell", areaCell, __LINE__); 
     io.read<real1d>("triangleQuality", triangleQuality, __LINE__); 
     io.read<real1d>("triangleAngleQuality", triangleAngleQuality, __LINE__); 
     io.read<real1d>("cellQuality", cellQuality, __LINE__); 
     io.read<real1d>("gridSpacing", gridSpacing, __LINE__); 
     io.read<real1d>("meshDensity", meshDensity, __LINE__); 

     io.read<int2d>("edgesOnEdge", edgesOnEdge, __LINE__); 
     io.read<int2d>("cellsOnEdge", cellsOnEdge, __LINE__); 
     io.read<int2d>("verticesOnEdge", verticesOnEdge, __LINE__); 
     io.read<int2d>("cellsOnVertex", cellsOnVertex, __LINE__); 
     io.read<int2d>("edgesOnVertex", edgesOnVertex, __LINE__); 
     io.read<int2d>("cellsOnCell", cellsOnCell, __LINE__); 
     io.read<int2d>("edgesOnCell", edgesOnCell, __LINE__); 
     io.read<int2d>("verticesOnCell", verticesOnCell, __LINE__); 
     io.read<int1d>("nEdgesOnEdge", nEdgesOnEdge, __LINE__); 
     io.read<int1d>("nEdgesOnCell", nEdgesOnCell, __LINE__); 

     io.read<int1d>("indexToEdgeID", indexToEdgeID, __LINE__); 
     io.read<int1d>("indexToCellID", indexToCellID, __LINE__); 
     io.read<int1d>("indexToVertexID", indexToVertexID, __LINE__); 

     io.read<int1d>("boundaryVertex", boundaryVertex, __LINE__); 
     io.read<int1d>("obtuseTriangle", obtuseTriangle, __LINE__); 

     io.read<real1d>("fEdge", fEdge, __LINE__); 
     io.read<real1d>("fCell", fCell, __LINE__); 
     io.read<real1d>("fVertex", fVertex, __LINE__); 

     io.read<real1d>("bottomDepth", bottomDepth, __LINE__); 
     io.close();

     std::cout << "done reading mesh" << std::endl;

     fixIndices();
     computeEdgeSign();
     
  }

  template<int N>
  void copy_to_device(Mesh<N> &mesh) {

     mesh.nCells = nCells;
     mesh.nEdges = nEdges;
     mesh.nVertices = nVertices;
     mesh.nVertLevels = nVertLevels;
     mesh.maxEdges = maxEdges;
     mesh.maxEdges2 = maxEdges2;
     mesh.vertexDegree = vertexDegree;

     mesh.xCell = xCell.createDeviceCopy();
     mesh.yCell = yCell.createDeviceCopy();
     mesh.zCell = zCell.createDeviceCopy();
     mesh.lonCell = lonCell.createDeviceCopy();
     mesh.latCell = latCell.createDeviceCopy();

     mesh.xEdge = xEdge.createDeviceCopy();
     mesh.yEdge = yEdge.createDeviceCopy();
     mesh.zEdge = zEdge.createDeviceCopy();
     mesh.lonEdge = lonEdge.createDeviceCopy();
     mesh.latEdge = latEdge.createDeviceCopy();

     mesh.xVertex = xVertex.createDeviceCopy();
     mesh.yVertex = yVertex.createDeviceCopy();
     mesh.zVertex = zVertex.createDeviceCopy();
     mesh.lonVertex = lonVertex.createDeviceCopy();
     mesh.latVertex = latVertex.createDeviceCopy();

     mesh.weightsOnEdge = weightsOnEdge.createDeviceCopy();

     mesh.angleEdge = angleEdge.createDeviceCopy();
     mesh.dcEdge = dcEdge.createDeviceCopy();
     mesh.dvEdge = dvEdge.createDeviceCopy();

     mesh.kiteAreasOnVertex = kiteAreasOnVertex.createDeviceCopy();
     mesh.areaTriangle = areaTriangle.createDeviceCopy();
     mesh.areaCell = areaCell.createDeviceCopy();
     mesh.triangleQuality = triangleQuality.createDeviceCopy();
     mesh.triangleAngleQuality = triangleAngleQuality.createDeviceCopy();
     mesh.cellQuality = cellQuality.createDeviceCopy();
     mesh.gridSpacing = gridSpacing.createDeviceCopy();
     mesh.meshDensity = meshDensity.createDeviceCopy();

     mesh.edgesOnEdge = edgesOnEdge.createDeviceCopy();
     mesh.cellsOnEdge = cellsOnEdge.createDeviceCopy();
     mesh.verticesOnEdge = verticesOnEdge.createDeviceCopy();
     mesh.cellsOnVertex = cellsOnVertex.createDeviceCopy();
     mesh.edgesOnVertex = edgesOnVertex.createDeviceCopy();
     mesh.cellsOnCell = cellsOnCell.createDeviceCopy();
     mesh.edgesOnCell = edgesOnCell.createDeviceCopy();
     mesh.verticesOnCell = verticesOnCell.createDeviceCopy();
     mesh.nEdgesOnEdge = nEdgesOnEdge.createDeviceCopy();
     mesh.nEdgesOnCell = nEdgesOnCell.createDeviceCopy();

     mesh.indexToEdgeID = indexToEdgeID.createDeviceCopy();
     mesh.indexToCellID = indexToCellID.createDeviceCopy();
     mesh.indexToVertexID = indexToVertexID.createDeviceCopy();

     mesh.boundaryVertex = boundaryVertex.createDeviceCopy();
     mesh.obtuseTriangle = obtuseTriangle.createDeviceCopy();

     mesh.fEdge = fEdge.createDeviceCopy();
     mesh.fCell = fCell.createDeviceCopy();
     mesh.fVertex = fVertex.createDeviceCopy();

     mesh.bottomDepth = bottomDepth.createDeviceCopy();

     mesh.edgeSignOnCell = edgeSignOnCell.createDeviceCopy();

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
