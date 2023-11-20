#include "const.h"
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
   int2dHost indexToEdgeID, indexToVertexID, indexToCellID;
   int2dHost nEdgesOnEdge, nEdgesOnCell;
   int2dHost boundaryVertex;
   int2dHost obtuseTriangle;
  
   real2dHost weightsOnEdge;
   real2dHost xCell, yCell, zCell;
   real2dHost lonCell, latCell;
   real2dHost xEdge, yEdge, zEdge;
   real2dHost lonEdge, latEdge;
   real2dHost xVertex, yVertex, zVertex;
   real2dHost lonVertex, latVertex;
   real2dHost dcEdge, dvEdge;
   real2dHost angleEdge;
   real2dHost cellQuality, triangleQuality, triangleAngleQuality;
   real2dHost areaCell, areaTriangle;
   real2dHost gridSpacing, meshDensity;
   real2dHost kiteAreasOnVertex;

   void read(const char* mesh_file){
     std::cout << mesh_file << "\n";

     ncwrap(ncmpi_open(MPI_COMM_WORLD, mesh_file, NC_NOWRITE, MPI_INFO_NULL, &nc_id), __LINE__);

     nCells = read_dim("nCells", __LINE__);
     nEdges = read_dim("nEdges", __LINE__);
     nVertices = read_dim("nVertices", __LINE__);
     maxEdges = read_dim("maxEdges", __LINE__);
     maxEdges2 = read_dim("maxEdges2", __LINE__);
     vertexDegree = read_dim("vertexDegree", __LINE__);

     xCell = read_double_var("xCell", 1, nCells, __LINE__); 
     yCell = read_double_var("yCell", 1, nCells, __LINE__); 
     zCell = read_double_var("zCell", 1, nCells, __LINE__); 
     lonCell = read_double_var("lonCell", 1, nCells, __LINE__); 
     latCell = read_double_var("latCell", 1, nCells, __LINE__); 

     xEdge = read_double_var("xEdge", 1, nEdges, __LINE__); 
     yEdge = read_double_var("yEdge", 1, nEdges, __LINE__); 
     zEdge = read_double_var("zEdge", 1, nEdges, __LINE__); 
     lonEdge = read_double_var("lonEdge", 1, nEdges, __LINE__); 
     latEdge = read_double_var("latEdge", 1, nEdges, __LINE__); 

     xVertex = read_double_var("xVertex", 1, nVertices, __LINE__); 
     yVertex = read_double_var("yVertex", 1, nVertices, __LINE__); 
     zVertex = read_double_var("zVertex", 1, nVertices, __LINE__); 
     lonVertex = read_double_var("lonVertex", 1, nVertices, __LINE__); 
     latVertex = read_double_var("latVertex", 1, nVertices, __LINE__); 

     weightsOnEdge = read_double_var("weightsOnEdge", nEdges, maxEdges2, __LINE__); 

     angleEdge = read_double_var("angleEdge", 1, nEdges, __LINE__); 
     dcEdge = read_double_var("dcEdge", 1, nEdges, __LINE__); 
     dvEdge = read_double_var("dvEdge", 1, nEdges, __LINE__); 

     kiteAreasOnVertex = read_double_var("kiteAreasOnVertex", nVertices, vertexDegree, __LINE__); 
     areaTriangle = read_double_var("areaTriangle", 1, nVertices, __LINE__); 
     areaCell= read_double_var("areaCell", 1, nCells, __LINE__); 
     triangleQuality = read_double_var("triangleQuality", 1, nVertices, __LINE__); 
     triangleAngleQuality = read_double_var("triangleAngleQuality", 1, nVertices, __LINE__); 
     cellQuality = read_double_var("cellQuality", 1, nCells, __LINE__); 
     gridSpacing = read_double_var("gridSpacing", 1, nCells, __LINE__); 
     meshDensity = read_double_var("meshDensity", 1, nCells, __LINE__); 

     edgesOnEdge = read_int_var("edgesOnEdge", nEdges, maxEdges2, __LINE__); 
     cellsOnEdge = read_int_var("cellsOnEdge", nEdges, 2, __LINE__); 
     verticesOnEdge = read_int_var("verticesOnEdge", nEdges, 2, __LINE__); 
     cellsOnVertex = read_int_var("cellsOnVertex", nVertices, vertexDegree, __LINE__); 
     edgesOnVertex = read_int_var("edgesOnVertex", nVertices, vertexDegree, __LINE__); 
     cellsOnCell = read_int_var("cellsOnCell", nCells, maxEdges, __LINE__); 
     edgesOnCell = read_int_var("edgesOnCell", nCells, maxEdges, __LINE__); 
     verticesOnCell = read_int_var("verticesOnCell", nCells, maxEdges, __LINE__); 
     nEdgesOnEdge = read_int_var("nEdgesOnEdge", 1, nEdges, __LINE__); 
     nEdgesOnCell = read_int_var("nEdgesOnCell", 1, nCells, __LINE__); 

     indexToEdgeID = read_int_var("indexToEdgeID", 1, nEdges, __LINE__); 
     indexToCellID = read_int_var("indexToCellID", 1, nCells, __LINE__); 
     indexToVertexID = read_int_var("indexToVertexID", 1, nVertices, __LINE__); 

     boundaryVertex = read_int_var("boundaryVertex", 1, nVertices, __LINE__); 
     obtuseTriangle = read_int_var("obtuseTriangle", 1, nVertices, __LINE__); 


     
  }

private:   
  int nc_id;

  void ncwrap(int err, int line) {
    if (err != NC_NOERR) {
      printf("NetCDF Error at line: %d\n", line);
      printf("%s\n",ncmpi_strerror(err));
      exit(-1);  
    }
  }

  int read_dim(const char* dim_str, int line) {
    int dim_id;
    int dim;
    MPI_Offset dim_read;

    ncwrap(ncmpi_inq_dimid(nc_id, dim_str, &dim_id), line);
    ncwrap(ncmpi_inq_dimlen(nc_id, dim_id, &dim_read), line);
    dim = dim_read;
   
    std::cout << dim_str << ": " << dim << "\n";

    return dim;
  }

  real2dHost read_double_var(const char* var_str,  int dim1, int dim2, int line) {
    int var_id;
    int i, j, k;
    real2dHost var;

    ncwrap(ncmpi_inq_varid(nc_id, var_str, &var_id), line);
    double buff[dim1*dim2];
    ncwrap(ncmpi_get_var_double_all(nc_id, var_id, buff), line);
    var = real2dHost(var_str, dim1, dim2);

    std::cout << "\n";
    std::cout << var_str << "\n";

    k = 0;
    for (j=0; j<dim1; j++) {
      for (i=0; i<dim2; i++) {
        var(j,i) = buff[k];
        std::cout << var(j,i) << ", ";
	k++;
      }
      std::cout << "\n";
    }
    return var;
  }

  int2dHost read_int_var(const char* var_str,  int dim1, int dim2, int line) {
    int var_id;
    int i, j, k;
    int2dHost var;

    ncwrap(ncmpi_inq_varid(nc_id, var_str, &var_id), line);
    int buff[dim1*dim2];
    ncwrap(ncmpi_get_var_int_all(nc_id, var_id, buff), line);
    var = int2dHost(var_str, dim1, dim2);

    std::cout << "\n";
    std::cout << var_str << "\n";

    k = 0;
    for (j=0; j<dim1; j++) {
      for (i=0; i<dim2; i++) {
        var(j,i) = buff[k];
        std::cout << var(j,i) << ", ";
	k++;
      }
      std::cout << "\n";
    }
    return var;
  }

};	
