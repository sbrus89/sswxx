#include "const.h"

class IO{
public:

  int nc_id;

  int t_dimid;
  int cell_dimid;
  int edge_dimid;
  int layer_dimid;
  int layerThickness_varid;
  int normalVelocity_varid;
  int ssh_varid;

  int output_count;

  void open(const char* nc_file) {
    ncwrap(ncmpi_open(MPI_COMM_WORLD, nc_file, NC_NOWRITE, MPI_INFO_NULL, &nc_id), __LINE__);
    //ncwrap(nc_open(nc_file, NC_NOWRITE, &nc_id), __LINE__);
  }

  void close() {
    ncwrap(ncmpi_close(nc_id), __LINE__);
    //ncwrap(nc_close(nc_id), __LINE__);
  }

  void ncwrap(int err, int line) {
    if (err != NC_NOERR) {
      printf("NetCDF Error at line: %d\n", line);
      printf("%s\n",ncmpi_strerror(err));
      //printf("%s\n",nc_strerror(err));
      exit(-1);
    }
  }

  int read_dim(const char* dim_str, int line) {
    int dim_id;
    int dim;
    MPI_Offset dim_read;
    //size_t dim_read;

    ncwrap(ncmpi_inq_dimid(nc_id, dim_str, &dim_id), line);
    ncwrap(ncmpi_inq_dimlen(nc_id, dim_id, &dim_read), line);
    //ncwrap(nc_inq_dimid(nc_id, dim_str, &dim_id), line);
    //ncwrap(nc_inq_dimlen(nc_id, dim_id, &dim_read), line);
    dim = dim_read;

    std::cout << dim_str << ": " << dim << std::endl;

    return dim;
  }

  template <class T> 
  void read(const char* var_str, T &var, int line) {
    int var_id;
    int type;
    int ndims;
    //int dimids[NC_MAX_VAR_DIMS];
    int dimids[10];
    int dim_start;
    int i, j;
    int n=1;
    MPI_Offset dim_read;
    //size_t dim_read;

    std::cout << std::endl;
    std::cout << var_str << std::endl; 

    ncwrap(ncmpi_inq_varid(nc_id, var_str, &var_id), line);
    ncwrap(ncmpi_inq_vartype(nc_id, var_id, &type), line);
    ncwrap(ncmpi_inq_varndims(nc_id, var_id, &ndims), line);
    ncwrap(ncmpi_inq_vardimid(nc_id, var_id, dimids), line);

    //ncwrap(nc_inq_varid(nc_id, var_str, &var_id), line);
    //ncwrap(nc_inq_vartype(nc_id, var_id, &type), line);
    //ncwrap(nc_inq_varndims(nc_id, var_id, &ndims), line);
    //ncwrap(nc_inq_vardimid(nc_id, var_id, dimids), line);

    std::cout << "# dims: " << ndims << std::endl;

    if (ndims == 3) {
      dim_start = 1;
    }
    else { 
      dim_start = 0;
    }
    std::vector<int> dims(ndims-dim_start);

    j = 0;
    for (i=dim_start; i<ndims; i++) {
      ncwrap(ncmpi_inq_dimlen(nc_id, dimids[i], &dim_read), line);
      //ncwrap(nc_inq_dimlen(nc_id, dimids[i], &dim_read), line);
      dims[j] = dim_read;
      std::cout <<  dims[j] << std::endl;
      j = j + 1;
    }
    yakl::Dims yakl_dims(dims);
    var = T(var_str, yakl_dims);

    for (i=0; i<dims.size(); i++) {
      n = n * dims[i];
    }

    if (type == NC_DOUBLE) { 
      std::cout << "double" << std::endl;
      double buff[n];
      ncwrap(ncmpi_get_var_double_all(nc_id, var_id, buff), line);
      //ncwrap(nc_get_var_double(nc_id, var_id, buff), line);
      for (i=0; i<n; i++) {
        var.data()[i] = buff[i];
      }
    }
    else if (type == NC_INT) {
      std::cout << "int" << std::endl;
      int buff[n];
      ncwrap(ncmpi_get_var_int_all(nc_id, var_id, buff), line);
      //ncwrap(nc_get_var_int(nc_id, var_id, buff), line);
      for (i=0; i<n; i++) {
        var.data()[i] = buff[i];
      }
    }
    else {
      printf("Error at line: %d\n", line);
      printf("Unsupported type\n");
      exit(-1);
    }

    //print_array(var);
  }

  template <class T> 
  void create(const char* nc_file, T &mesh) {

    int dimids[3];

    ncwrap(ncmpi_create(MPI_COMM_WORLD, nc_file, NC_NOWRITE, MPI_INFO_NULL, &nc_id), __LINE__);

    ncwrap(ncmpi_def_dim(nc_id, "Time", (MPI_Offset) NC_UNLIMITED, &t_dimid), __LINE__);
    ncwrap(ncmpi_def_dim(nc_id, "nCells", (MPI_Offset) mesh.nCells, &cell_dimid), __LINE__);
    ncwrap(ncmpi_def_dim(nc_id, "nEdges", (MPI_Offset) mesh.nEdges, &edge_dimid), __LINE__);
    ncwrap(ncmpi_def_dim(nc_id, "nVertLayers", (MPI_Offset) mesh.nVertLevels, &layer_dimid), __LINE__);

    //ncwrap(nc_create(nc_file, NC_NOWRITE, &nc_id), __LINE__);

    size_t nCells, nEdges, nVertLevels;
    nCells = mesh.nCells;
    nEdges = mesh.nEdges;
    nVertLevels = mesh.nVertLevels;

    //ncwrap(nc_def_dim(nc_id, "Time", NC_UNLIMITED, &t_dimid), __LINE__);
    //ncwrap(nc_def_dim(nc_id, "nCells",  nCells, &cell_dimid), __LINE__);
    //ncwrap(nc_def_dim(nc_id, "nEdges",  nEdges, &edge_dimid), __LINE__);
    //ncwrap(nc_def_dim(nc_id, "nVertLayers", nVertLevels, &layer_dimid), __LINE__);

    dimids[0] = t_dimid;
    dimids[1] = cell_dimid;
    dimids[2] = layer_dimid;

    ncwrap(ncmpi_def_var(nc_id, "layerThickness", NC_DOUBLE, 3, dimids, &layerThickness_varid), __LINE__);
    ncwrap(ncmpi_def_var(nc_id, "ssh", NC_DOUBLE, 2, dimids, &ssh_varid), __LINE__);
    //ncwrap(nc_def_var(nc_id, "layerThickness", NC_DOUBLE, 3, dimids, &layerThickness_varid), __LINE__);
    //ncwrap(nc_def_var(nc_id, "ssh", NC_DOUBLE, 2, dimids, &ssh_varid), __LINE__);
    dimids[1] = edge_dimid;
    ncwrap(ncmpi_def_var(nc_id, "normalVelocity", NC_DOUBLE, 3, dimids, &normalVelocity_varid), __LINE__);
    //ncwrap(nc_def_var(nc_id, "normalVelocity", NC_DOUBLE, 3, dimids, &normalVelocity_varid), __LINE__);

    ncwrap(ncmpi_enddef(nc_id), __LINE__);
    //ncwrap(nc_enddef(nc_id), __LINE__);

    output_count = 0;

  }

  template <class T> 
  void write(T &state) {

    MPI_Offset st3[3], ct3[3];
    //std::vector<size_t> st3(3), ct3(3);

    st3[0] = output_count;
    st3[1] = 0;
    st3[2] = 0;

    ct3[0] = 1;
    ct3[1] = state.nCells;
    ct3[2] = state.nVertLevels;

    ncwrap(ncmpi_put_vara_double_all(nc_id,  layerThickness_varid, st3, ct3, state.layerThickness.data()) , __LINE__);
    ncwrap(ncmpi_put_vara_double_all(nc_id,  ssh_varid, st3, ct3, state.ssh.data()) , __LINE__);
    //ncwrap(nc_put_vara_double(nc_id,  layerThickness_varid, st3.data(), ct3.data(), state.layerThickness.data()) , __LINE__);
    //ncwrap(nc_put_vara_double(nc_id,  ssh_varid, st3.data(), ct3.data(), state.ssh.data()) , __LINE__);
    ct3[1] = state.nEdges;
    ncwrap(ncmpi_put_vara_double_all(nc_id,  normalVelocity_varid, st3, ct3, state.normalVelocity.data()) , __LINE__);
    //ncwrap(nc_put_vara_double(nc_id,  normalVelocity_varid, st3.data(), ct3.data(), state.normalVelocity.data()) , __LINE__);

    output_count = output_count + 1;
  }

  void print_array(int1dHost var) {
    int j;
    
    for (j=0; j<var.dimension[0]; j++) {
        std::cout << var(j) << "\n";
    }
  }

  void print_array(real1dHost var) {
    int j;
    
    for (j=0; j<var.dimension[0]; j++) {
        std::cout << var(j) << "\n";
    }
  }

  void print_array(int2dHost var) {
    int i, j;
    
    for (j=0; j<var.dimension[0]; j++) {
      for (i=0; i<var.dimension[1]; i++) {
        std::cout << var(j,i) << ", ";
      }
      std::cout << "\n";
    }
  }

  void print_array(real2dHost var) {
    int i, j;
    std::cout << var.dimension[0] << "\n";
    std::cout << var.dimension[1] << "\n";
    for (j=0; j<var.dimension[0]; j++) {
      for (i=0; i<var.dimension[1]; i++) {
        std::cout << var(j,i) << ", ";
      }
      std::cout << "\n";
    }
  }
};
