#include "const.h"

class IO{
public:

  int nc_id;

  void open(const char* nc_file) {
    ncwrap(ncmpi_open(MPI_COMM_WORLD, nc_file, NC_NOWRITE, MPI_INFO_NULL, &nc_id), __LINE__);
  }

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

  template <class T> 
  T read(const char* var_str, int line) {
    int var_id;
    int type;
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    int i;
    int n=1;
    T var;
    MPI_Offset dim_read;

    std::cout << "\n";
    std::cout << var_str << "\n"; 

    ncwrap(ncmpi_inq_varid(nc_id, var_str, &var_id), line);
    ncwrap(ncmpi_inq_vartype(nc_id, var_id, &type), line);
    ncwrap(ncmpi_inq_varndims(nc_id, var_id, &ndims), line);
    ncwrap(ncmpi_inq_vardimid(nc_id, var_id, dimids), line);

    std::cout << "# dims: " << ndims << "\n";

    std::vector<int> dims(ndims);
    for (i=0; i<ndims; i++) {
      ncwrap(ncmpi_inq_dimlen(nc_id, dimids[i], &dim_read), line);
      dims[i] = dim_read;
      std::cout <<  dims[i] << "\n";
    }


    yakl::Dims yakl_dims(dims);
    var = T(var_str, yakl_dims);

    for (i=0; i<dims.size(); i++) {
      n = n * dims[i];
    }

    if (type == NC_DOUBLE) { 
      std::cout << "double\n";
      double buff[n];
      ncwrap(ncmpi_get_var_double_all(nc_id, var_id, buff), line);
      for (i=0; i<n; i++) {
        var.data()[i] = buff[i];
      }
    }
    else if (type == NC_INT) {
      std::cout << "int\n";
      int buff[n];
      ncwrap(ncmpi_get_var_int_all(nc_id, var_id, buff), line);
      for (i=0; i<n; i++) {
        var.data()[i] = buff[i];
      }
    }
    else {
      printf("Error at line: %d\n", line);
      printf("Unsupported type\n");
      exit(-1);
    }

    print_array(var);

    return var;

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
    
    for (j=0; j<var.dimension[0]; j++) {
      for (i=0; i<var.dimension[1]; i++) {
        std::cout << var(j,i) << ", ";
      }
      std::cout << "\n";
    }
  }
};
