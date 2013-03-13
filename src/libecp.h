/* Copyright (c) 2012, Christoph Reimann */

#ifndef LIBECP_H
#define LIBECP_H 1

typedef struct _libECPHandle libECPHandle; 

typedef void (*ECPCallback)(int A, int B, int C, 
			    int sa, int sb, 
			    int la, int lb, 
			    int shifta, int shiftb,
			    double *I, void *p);

/* create libECP handle and initialize type1 and type2 integrations */
libECPHandle * libECP_init(int nrAtoms, double *geometry, 
			   int *shellsECP, int *KECP, 
			   int *lECP, double *nECP, double *dECP, double *aECP, 
			   int *shellsBS, int *lBS, int *KBS, 
			   double *dBS, double *aBS, 
			   int n, int lmax, int *shellOrdering, 
			   int largeGridOrder, double tolerance,  double accuracy);

/* check return value:
   0 - successful integration
   1 - error during type1 integration
   2 - error during type2 integration */
int calculateECPIntegrals(libECPHandle *h, ECPCallback cb, void *args);

void libECP_free(libECPHandle *h);

#endif
