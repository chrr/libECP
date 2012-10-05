
#ifndef TYPE1_H
#define TYPE1_H 1

#include "gc_integrators.h"
#include "bessel.h"
#include "ecp.h"

typedef struct {
  /* number of Gaussian primitives per shell */
  int *contraction;
  /* coefficients d and exponents a of primitive Gaussians */
  double *d;
  double *a;
  /* factorials */
  double *fac;
  double *dfac;
  /* shift: order of basis function derivative w.r.t. nuclear coordinates */
  int shift;
  int *ijk;
  int *ijkIndex;
  int ijkDim;
  /* transformation matrices */
  double *poly2sph;
  int tmDim;

  BesselFunctionTable *bessel;
  int maxLambda;
  int largeOrder;
  GCIntegrationTable *largeGrid;
  int smallOrder;
  GCIntegrationTable *smallGrid;
  /* general accuracy parameter - "numerical zero" */
  double accuracy;
  double lnAccuracy;
  double tolerance;
} Type1;


/* allocate memory for Type2 structure and set default dimensions */
Type1 * Type1_new(int nrAtoms, double *geometry, ECP **U, 
		  int *shells, int *am, int *contraction, double *d, double *a,
		  int derivative, const double tolerance);

void Type1_init(Type1 *t1, const int maxLambda, double *fac, double *dfac, int *ijk, int *ijkIndex, const int ijkDim,
		double *poly2sph, const int tmDim, BesselFunctionTable *bessel, int largeOrder, GCIntegrationTable *largeGrid, 
		const double accuracy, const double lnAccuracy);

void Type1_free(Type1 *t1);

doubleArray calcChi(Type1 *t1, double *U_L, ECP *U,
		    double *rAC, const double dAC, double *rBC, const double dBC, 
		    const int la, const int shella, const int offseta, const int shifta, doubleArray uspA,  
		    const int lb, const int shellb, const int offsetb, const int shiftb, doubleArray uspB);

#endif
