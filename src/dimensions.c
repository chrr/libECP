/* Copyright (c) 2012, Christoph Reimann */

#include <stdlib.h>

#include "dimensions.h"


/* get exponents nx,ny,nz for cartesian gaussian shells with
   angular momentum up to am
   default: standard LIBINT orbital ordering, i.e.
            p : px , py , pz
	    d : dxx , dxy , dxz , dyy , dyz , dzz
	    f : fxxx , fxxy , fxxz , fxyy , fxyz , fxzz , fyyy , fyyz , fyzz , fzzz
	    etc. 
   parameter: max. angular momentum 
              special case am = -1: reset static variables */
int *
cartesianShellOrder(const int am) {
  int i,j, l, index, c, nx;
  /* x,y,z exponents */
  int *indices = calloc(3*C_DIM(am), sizeof(int));
  
  for (l=0; l<=am; l++) {
    c = 0;
    /* determine set of exponents c for shell with angular momentum l */
    for(i=0; i<=l; i++) {
      nx = l-i; /* nx */
      for(j=0; j<=i; j++) {
	index = CIJK_INDEX(l,c);
	indices[index]   = nx; 
	indices[index+1] = i-j;  /* ny */
	indices[index+2] = j;    /* nz */
	c++;
      }
    }
  }    
  return indices;
}


int *
cartesianShellOrderIndex(const int am, int *ijk) {
  int i, j, k, l, index, c;
  const int dim = am+1;
  int *indices = calloc(dim*dim*dim, sizeof(int));

  for (l=0; l<=am; l++) {
    for(c=0; c<IJK_DIM(l); c++) {
      index = CIJK_INDEX(l,c);
      i = ijk[index];
      j = ijk[index+1];
      k = ijk[index+2];
      indices[i*dim*dim+j*dim+k] = index/3;
    }
  }    
  return indices;
}

