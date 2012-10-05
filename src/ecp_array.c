
#include "ecp_array.h"

#include <stdarg.h>
#include <assert.h>


doubleArray
doubleArray_new(size_t dim, ...) {
  va_list args;
  int i;
  doubleArray array;
  size_t dataSize = 1;

  array.dim = calloc(dim+1, sizeof(size_t));
  array.dim[0] = dim;
  
  /* get the row dimensions */
  va_start(args, dim);
  for(i=1; i<=dim; i++) {
    array.dim[i] = va_arg(args, size_t);
    assert(array.dim[i] > 0);
    dataSize *= array.dim[i];
  }

  array.array = calloc(dataSize, sizeof(double));
  for(i=0; i<dataSize; i++)
    array.array[i] = 0.0;
  return array;
}


void
doubleArray_free(doubleArray arr) {
  free(arr.dim);
  free(arr.array);
  arr.dim = NULL;
  arr.array = NULL;
}

