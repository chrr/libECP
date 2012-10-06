/* Copyright (c) 2012, Christoph Reimann */

#ifndef LIBECP_H
#define LIBECP_H 1

typedef struct _libECPHandle libECPHandle; 

typedef void (*ECPCallback)(int A, int sa, int la, int shifta, 
			    int B, int sb, int lb, int shiftb,
			    int C, double *I, void *p);

libECPHandle * libECP_init(int nrAtoms, double *geometry, 
			   int *shellsECP, int *contractionECP, int *lECP, double *nECP, double *dECP, double *aECP, 
			   int *shells, int *am, int *contraction, double *d, double *a, 
			   int derivative, int *shellOrdering, const int lmax);

/* check return value:
   0 - successful integration
   1 - error during type1 integration
   2 - error during type2 integration */
int calculateECPIntegrals(libECPHandle *h, ECPCallback cb, void *args);

void libECP_free(libECPHandle *h);

#endif
