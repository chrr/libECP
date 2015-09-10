/* Copyright (c) 2013, Christoph Reimann */

#ifndef GET_INTEGRALS_H
#define GET_INTEGRALS_H 1


int getIntegrals (int nrAtoms, double *geometry, 
		  int *shellsECP, int *KECP, 
		  int *lECP, double *nECP, double *dECP, double *aECP, 
		  int *shellsBS, int *lBS, int *KBS, 
		  double *dBS, double *aBS, 
		  int largeGridOrder, double tolerance,  double accuracy, 
		  int rowdim, double *I);

#endif /* GET_INTEGRALS_H */
