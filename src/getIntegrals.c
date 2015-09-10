
#include <stdlib.h>
#include <stdio.h>

#include "libecp.h"
#include "getIntegrals.h"
#include "dimensions.h"

typedef struct {
  /* integral matrix is in fact symmetric, so we waste some memory here */
  double *I;
  int dimI;
  int *aoDim;
  int dim;
  /* for use in gradient calculations
  double **P;
  double *gradient;
  */
} libECP_params;


void
libECP_callback0(int A, int sa, int la, int shifta,
		 int B, int sb, int lb, int shiftb,
		 int C, double *I, void *p) {
  int i, j;
  libECP_params *P = (libECP_params *)p;
  double *integrals = P->I;
  const int dimI = P->dimI;
  const int *aodim = P->aoDim;
  const int dim = P->dim;
  const int lstartA = aodim[A*dim+sa];
  const int lstartB = aodim[B*dim+sb];

  /* copy matrix elements */
  for(i=0; i<IJK_DIM(la); i++) {
    for(j=0; j<IJK_DIM(lb); j++) {
      if(lstartA+i>lstartB+j)
	continue;
      integrals[(lstartA+i)*dimI+(lstartB+j)] += I[i*IJK_DIM(lb)+j];
    }
  }
}

int getIntegrals (int nrAtoms, double *geometry, 
		  int *shellsECP, int *KECP, 
		  int *lECP, double *nECP, double *dECP, double *aECP, 
		  int *shellsBS, int *lBS, int *KBS, 
		  double *dBS, double *aBS, 
		  /* 1024, 1.0E-12, 1.0E-14 */
		  int largeGridOrder, double tolerance,  double accuracy, 
		  int rowdim, double *I) {
  libECP_params *P = malloc(sizeof(libECP_params));
  libECPHandle *h;
  int *aoDim, dim; 
  int maxShells = 0;
  int i, j, lstart = 0;
  int idx = 0;
  /* initialize aoDim */
  for(i=0; i<nrAtoms; i++) {
    if (shellsBS[i] > maxShells)
      maxShells = shellsBS[i];
  }
  dim = maxShells;
  aoDim = calloc(nrAtoms*dim, sizeof(int));
  for (i=0; i<nrAtoms; i++) {
    for (j=0; j<shellsBS[i]; j++) {
      aoDim[i*dim+j] = lstart;
      lstart += IJK_DIM(lBS[idx]);
      idx++;
    }
  }
  P->I = I;
  P->dimI = rowdim;
  P->aoDim = aoDim;
  P->dim = dim;

  h = libECP_init(nrAtoms, geometry, shellsECP, lECP, KECP, nECP, dECP, aECP,
		  shellsBS, lBS, KBS, dBS, aBS,
		  0, /* order for derivative calculations; 0 == normal integral run */
		  -1, NULL, /* use default (== libint) shell ordering */
		  /* parameters determining the accuracy of the numerical integration */
		  largeGridOrder, tolerance, accuracy); 
  if(NULL==h) {
    printf("error initializing libECP\n");
    return 1;
  }
  calculateECPIntegrals(h, libECP_callback0, P);
  libECP_free(h);

  free(aoDim);
  free(P);

  return 0;  
}
