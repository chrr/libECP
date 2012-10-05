
#ifndef UTIL_H
#define UTIL_H 1

#include "ecp_array.h"
#include "gc_integrators.h"

typedef enum {
  ECP_FALSE = 0,
  ECP_TRUE = 1
} ECPBool;

/* calculate n! ineffectively */
double * factorials(const int n);

/* calculate n!! ineffectively (even and odd numbers n) */
double * double_factorials(const int n);

/* calculate the binomial coefficient n over k */
double n_over_k(const int n, const int k, const double *fac);

void sphericalCoordinates (double *cartesian, double *spherical);
double distance(double *A, double *B);
void distanceVector(double *AB, double *A, double *B);

/* calculate atomic radius for primitive gaussian function */
double gtoRadius (const double c, const double zeta, const int l, const double cutoff);
double shellRadius(int depth, int am, double *d, double *a, const double zero);

/* screening by potential cutoff
   returns: index (elements U[0..index] are to be retained)
   -1: potential is screened completely */
int potentialScreening(GCIntegrationTable *grid, double *U, double zero);

/* evaluate unitary sphere polynomial terms x^i y^j z^k */
doubleArray unitarySpherePolynomials(const int l, double *r);

doubleArray calcPolynomials(doubleArray gamma, const int *ijk, const int *ijkIndex, const int ijkDim, 
			    double *fac, const double tol, const double N, 
			    const int la, doubleArray uspA, const int lb, doubleArray uspB);

#endif
