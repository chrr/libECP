/* Copyright (c) 2012, Christoph Reimann */

#ifndef ECP_ARRAY_H
#define ECP_ARRAY_H 1

#include <stdlib.h>

typedef struct {
  size_t *dim;
  double *array; 
} doubleArray;

typedef struct {
  size_t *dim;
  int *array; 
} intArray;

doubleArray doubleArray_new(size_t dim, ...);
void doubleArray_free(doubleArray arr);

#endif
