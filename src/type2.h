
#ifndef TYPE2_H
#define TYPE2_H 1

#include "ecp.h"
#include "bessel.h"
#include "gc_integrators.h"
#include "ecp_array.h"
#include "util.h"

typedef enum {
  TYPE2_OK = 0, 
  BESSEL_NOT_CONVERGED = 1, 
} ErrType2;

typedef struct _Type2 Type2;

struct _Type2 {
  int nrAtoms;
  /* positions of the nuclei */ 
  double *geometry;
  /* basis set */
  /* number of shells per atom */
  int *shells;
  /* total number of shells := sum(shells) */
  int nrShells;
  /* maximum number of shells per atom */
  int maxShells;
  /* angular momentum of each shell */
  int *am;
  /* number of Gaussian primitives per shell */
  int *contraction;
  /* coefficients d and exponents a of primitive Gaussians */
  double *d;
  double *a;
  /* effective core potentials */
  ECP *U;
  int maxL;
  int maxLBS;
  int maxLambda;
  int maxAlpha;
  /* shift: order of basis function derivative w.r.t. nuclear coordinates */
  int shift;
  int *ijk;
  int *ijkIndex;
  int ijkDim;
  /* transformation matrices */
  double *cart2sph;
  double *poly2sph;
  int tmDim;

  /* factorials */
  double *fac;
  double *dfac;

  BesselFunctionTable *bessel;
  /* number of abscissas for bessel function tabulation */
  int N;
  /* upper limit for lambda argument */
  int lmax;
  /* cutoff for series expansion */
  int cutoff;
  /* grid orders for radial integration */
  int smallOrder;
  GCIntegrationTable *smallGrid;
  int largeOrder;
  GCIntegrationTable *largeGrid;
  /* tolerance used in radial integration */
  double tolerance;
  /* general accuracy parameter - "numerical zero" */
  double accuracy;
  double lnAccuracy;

  ECPBool freeMemory;
};

/* allocate memory for Type2 structure and set default dimensions */
Type2 * Type2_new(const int nrAtoms, double *geometry, ECP **U, 
		  int *shells, int *am, 
		  int *contraction, double *d, double *a, 
		  const int derivative, const double tol);

/* initialize and allocate memory for Type2 pointer fields */
void Type2_init(Type2 *t2, double *fac, double *dfac, int *ijk, int *ijkIndex, const int ijkDim, 
		double *cart2sph, double *poly2sph, const int tmDim, GCIntegrationTable *largeGrid, 
		BesselFunctionTable *bessel);
ErrType2 Type2_initDefault(Type2 *t2);

void Type2_free(Type2 *t2);


typedef struct {
  int start;
  int end;
  ECPBool skipShell;
} ScreenedGrid;


doubleArray calcF_FM06(Type2 *t2, double *r_C, ScreenedGrid **screening);
doubleArray tabECP_FM06(ECP *U, GCIntegrationTable *grid, const int lmax, ECPBool useScreening, const double tol);
doubleArray calcGamma(Type2 *t2, doubleArray F, doubleArray UTab, ECP *U, const double dAC, const double dBC, 
		      const int la, const int shella, const int shellOffseta, doubleArray omegaA, doubleArray uspA, const int shifta,  
		      const int lb, const int shellb, const int shellOffsetb, doubleArray omegaB, doubleArray uspB, const int shiftb);

#endif 
