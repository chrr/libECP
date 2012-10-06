/* Copyright (c) 2012, Christoph Reimann */

#ifndef ANGULAR_INTEGRAL_H
#define ANGULAR_INTEGRAL_H 1

#include "type2.h"

/* tabulate the angular integral Omega(l,m,lambda,mu,i,j,k) 
   which is independent of basis set and ECP */
doubleArray tabOmega(Type2 *t2); 

doubleArray evalAngularIntegrals(Type2 *t2, doubleArray Omega, int lmax_A, double *r_A);

#endif 
