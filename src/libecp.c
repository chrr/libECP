/* Copyright (c) 2012, Christoph Reimann */

#include <math.h>

#include "ecp.h"
#include "ecp_array.h"
#include "type1.h"
#include "type2.h"
#include "util.h"
#include "angular_integrals.h"
#include "libecp.h"
#include "dimensions.h"
#include "transformations.h"
#include "gc_integrators.h"


struct _libECPHandle {
  int nrAtoms;
  double *geometry;
  ECP **U;
  int *shells, *am, *contraction;
  double *d, *a;
  int derivative;
  int maxLECP;
  int maxLBS;
  int maxAlpha;
  int *maxLAtom;
  Type1 *t1;
  Type2 *t2;
  BesselFunctionTable *bessel;
  /* Gauss-Chebyshev grid */
  GCIntegrationTable *largeGrid;
  /* factorials */
  double *fac;
  double *dfac;
  int *ijk;
  int *ijkIndex;
  int ijkDim;
  int freeIJK;
  /* transformation matrices */
  int tmDim;
  double *cart2sph;
  double *poly2sph;
};


libECPHandle * 
libECP_init(int nrAtoms, double *geometry, 
	    int *shellsECP, int *contractionECP, int *lECP, double *nECP, double *dECP, double *aECP, 
	    int *shells, int *am, int *contraction, double *d, double *a, 
	    int derivative, int *shellOrdering, const int lmax) {
  libECPHandle *h = malloc(sizeof(libECPHandle));
  int i, j, k, L, N, index, maxLECP, maxLBS, maxAlpha, maxLambda;
  ECP **U = (ECP**)calloc(nrAtoms, sizeof(ECP*));
  int largeGridOrder = 1024;
  /* tolerance in quadratures for radial integrals */
  const double tolerance = 1.0E-14;
  /* general accuracy parameters - "numerical zero" */
  const double accuracy = 1.0E-14;
  const double lnAccuracy = log(accuracy)-2;

  /* initialize libECPHandle */
  h->nrAtoms = nrAtoms;
  h->geometry = geometry;
  h->U = U;
  h->shells = shells;
  h->am = am;
  h->contraction = contraction;
  h->d = d;
  h->a = a;
  h->derivative = derivative;
  h->maxLAtom = NULL;
  h->t1 = NULL;
  h->t2 = NULL;
  h->bessel = NULL;
  h->largeGrid = NULL;
  h->fac = NULL;
  h->dfac = NULL;
  h->ijk = NULL;
  h->ijkIndex = NULL;
  h->cart2sph = NULL;
  h->poly2sph = NULL;

  /* initialize ECP and determine maxLECP */
  index = 0;  
  maxLECP = 0;
  for(i=0; i<nrAtoms; i++) {
    U[i] = NULL;
    if(0<shellsECP[i]) {
      L = N = 0;
      for(j=0; j<shellsECP[i]; j++) {
	if (L<lECP[index])
	  L = lECP[index];
	N += contractionECP[index];
	index++;
      }
      U[i] = newECP(L, N);
      if (L>maxLECP)
	maxLECP = L;
    }
  }
  /* set ECP parameters */
  L = N = 0;
  for(i=0; i<nrAtoms; i++) {
    if(NULL!=U[i]) {
      index = 0;
      for(j=0; j<shellsECP[i]; j++) {
	for(k=0; k<contractionECP[L]; k++) { 
	  U[i]->U_l[index].l = lECP[L];
	  U[i]->U_l[index].n = nECP[N];
	  U[i]->U_l[index].d = dECP[N];
	  U[i]->U_l[index].a = aECP[N];
	  N++; index++;
	}
	L++;
      }
    }
  }

  /* determine maxLBS and maxLAtom */
  h->maxLAtom = calloc(nrAtoms, sizeof(int));
  index = 0;
  maxLBS = 0;
  for(i=0; i<nrAtoms; i++) {
    for(j=0; j<shells[i]; j++) {
      if (h->maxLAtom[i] < am[index])
	h->maxLAtom[i] = am[index];
      if (maxLBS < am[index])
	maxLBS = am[index];
      index++;
    }
  }

  maxAlpha = maxLBS + derivative;
  /* limits for lambda due to angular integration */ 
  maxLambda = maxLECP - 1 + maxAlpha;
  /* limit set by call to TM_DIM */
  h->fac = factorials(2*(maxLambda+maxAlpha)+1);
  /* limit set by call in poly2sph */
  h->dfac = double_factorials(2*(maxLambda+maxAlpha)+1);
  /* limits for ijk and ijkIndex set by transformation matrices */
  h->freeIJK = 0;
  if (NULL==shellOrdering) {
    h->freeIJK = 1;
    h->ijk = cartesianShellOrder(maxLambda+maxAlpha);
    h->ijkIndex = cartesianShellOrderIndex(maxLambda+maxAlpha, h->ijk);
    h->ijkDim = maxLambda+maxAlpha+1;
  }
  else {
    if(lmax<maxLambda+maxAlpha+1) {
      libECP_free(h);
      return NULL;
    }
    h->ijk = shellOrdering;
    h->ijkIndex = cartesianShellOrderIndex(lmax, h->ijk);
    h->ijkDim = lmax+1;
  }
  /* type2: maxLambda+maxAlpha
     type1: maxLBS+maxAlpha */
  h->tmDim = maxLambda+maxAlpha;
  h->cart2sph = TM_cart2sph(h->tmDim, h->ijk, h->fac);
  h->poly2sph = TM_poly2sph(h->cart2sph, h->tmDim, h->ijk, h->dfac);

  /* Gauss-Chebyshev grid */
  h->largeGrid = initGridGC_PSM92(largeGridOrder, tolerance);

  /* tabulate bessel function 
     upper limit for lambda argument: maxLECP+maxAlpha+6
     number of abscissas for bessel function tabulation: 16*100
     upper limit for series expansion: 200 */
  h->bessel = malloc(sizeof(BesselFunctionTable));
  if (0 != tabBesselFunction(h->bessel, maxLECP+maxAlpha+6, 16*100, 200, accuracy)) {
    /* error: no convergence */
    libECP_free(h);
    return NULL;
  }
  
  /* initialize parameters for type1 and type2 integrations */
  h->t1 = Type1_new(nrAtoms, geometry, U, shells, am, contraction, d, a, derivative, tolerance);
  Type1_init(h->t1, maxLambda, h->fac, h->dfac, h->ijk, h->ijkIndex, h->ijkDim, 
	     h->poly2sph, h->tmDim, h->bessel, largeGridOrder, h->largeGrid, accuracy, lnAccuracy);
  h->t2 = Type2_new(nrAtoms, geometry, U, shells, am, contraction, d, a, derivative, tolerance);
  Type2_init(h->t2, h->fac, h->dfac, h->ijk, h->ijkIndex, h->ijkDim, h->cart2sph, h->poly2sph, h->tmDim, 
	     h->largeGrid, h->bessel);
  /* if (TYPE2_OK != Type2_initDefault(h->t2)) { */
  /*   /\* error during type2 initialization *\/ */
  /*   libECP_free(h); */
  /*   return NULL; */
  /* } */

  return h;
}

typedef struct {
  int order;
  int N;
  /* in case of higher order derivatives,
     the hard-coded array limits have to be increased! */
  int shifta[4];
  int shiftb[4];
} derivativeShift;

int
calculateECPIntegrals(libECPHandle *h, ECPCallback cb, void *args) {
  Type1 *t1 = h->t1;
  Type2 *t2 = h->t2;
  const int nrAtoms = h->nrAtoms;
  double *geometry = h->geometry;
  const int d = h->derivative;
  int *atomMaxl = h->maxLAtom;
  ECP **U = h->U;
  int A, B, C;
  int la, lb, s1, s2;
  int s, shifta, shiftb;
  int gridStart, gridEnd;
  double *rA, *rB, *rC;
  double rAC[3], rBC[3], dAC, dBC;
  doubleArray omega, omegaA, omegaB, FTab, UTab, uspA, uspB;
  double *U_L = calloc(t1->smallGrid->order, sizeof(double));
  doubleArray chi, gamma;
  /* screening with respect to the basis set */
  ScreenedGrid *sg = NULL;
  /* indices for shells, orbitals and contractions of primitive Gaussians */
  int sOffset1, sOffset2, sIdx1, sIdx2;
  int pOffset1, pOffset2, pIdx1, pIdx2;
  /* derivatives: use translational invariance criterium to avoid 
     calculating derivatives of the potential itself
     e.g. first derivative of a primitive gaussian function g at atom A: 
       dg(ax,ay,az)/dAx = -ax g(ax-1,ay,ay) + 2 zeta g(ax+1,ay,az)
     i.e. derivatives involve gaussians shifted in angular momentum */			    
  const derivativeShift shifts[] = { { 0, 1, {0}, {0} }, 
				     { 1, 4, {+1, -1, 0, 0}, {0, 0, +1, -1} }, /* first derivative */ 
  };
  int result = 0;

  /* tabulate angular type2 integrals */
  omega = tabOmega(t2);

  for (C=0; C<nrAtoms; C++) {
    if(NULL == U[C]) {
      continue;
    }
    rC = geometry + C*3;

    /* tabulate ECP potential and F function values for fast type2 radial quadrature */
    t2->smallGrid->start = 0;
    t2->smallGrid->end = t2->smallGrid->order-1;
    UTab = tabECP_FM06(U[C], t2->smallGrid, t2->maxLambda, ECP_TRUE, t2->accuracy);
    FTab = calcF_FM06(t2, rC, &sg);

    /* calculate the local potential U_L at every grid point */
    for(s=0; s<t1->smallGrid->order; s++) 
      U_L[s] = evalECP(U[C], U[C]->L, t1->smallGrid->x[s]);
    t1->smallGrid->end = potentialScreening(t1->smallGrid, U_L, t1->accuracy);
    
    /* loop over atomic centers A,B */
    sIdx1 = pIdx1 = 0;
    sOffset1 = pOffset1 = 0;
    sIdx2 = pIdx2 = 0;
    sOffset2 = pOffset2 = 0;
    for(A=0; A<nrAtoms; A++) {
      rA = geometry + A*3;
      dAC = distance(rC, rA);
      distanceVector(rAC, rC, rA);
      uspA = unitarySpherePolynomials(atomMaxl[A]+d, rAC); 
      omegaA = evalAngularIntegrals(t2, omega, atomMaxl[A]+d, rAC);

      sOffset2 = sOffset1;
      pOffset2 = pOffset1;
      for(B=A; B<nrAtoms; B++) {
	rB = geometry + B*3;
	dBC = distance(rC, rB);
	distanceVector(rBC, rC, rB);
	uspB = unitarySpherePolynomials(atomMaxl[B]+d, rBC); 
	omegaB = evalAngularIntegrals(t2, omega, atomMaxl[B]+d, rBC);

	/* loop over shells */
	sIdx1 = sOffset1;
	pIdx1 = pOffset1;
	for(s1=0; s1<h->shells[A]; s1++) {
	  la = h->am[sIdx1];	  
	  sIdx2 = sOffset2;
	  pIdx2 = pOffset2;
	  for(s2=0; s2<h->shells[B]; s2++) {
	    lb = h->am[sIdx2];

	    if ((A==B && s2<s1) || 
		(1==d && A==C && B==C) || 
		(ECP_TRUE == sg[sIdx1].skipShell) ||
		(ECP_TRUE == sg[sIdx2].skipShell)) {
	      pIdx2 += h->contraction[sIdx2];
	      sIdx2++;
	      continue;
	    }

	    /* reset grid information - it might have been reduced by screening */   
	    gridStart = sg[sIdx1].start > sg[sIdx2].start ? sg[sIdx1].start : sg[sIdx2].start;
	    gridEnd = sg[sIdx1].end > sg[sIdx2].end ? sg[sIdx1].end : sg[sIdx2].end;
	    t1->smallGrid->start = gridStart;
	    t2->smallGrid->start = gridStart;   	    
	    t1->smallGrid->end = gridEnd;
	    t2->smallGrid->end = gridEnd;	    

	    for(s=0; s<shifts[d].N; s++) {
	      shifta = shifts[d].shifta[s];
	      shiftb = shifts[d].shiftb[s];
	      if ((-1==shifta && la<=0) || /* avoid negative angular momentum */ 
		  (-1==shiftb && lb<=0) ||
		  (C==A && 0!=shifta)   || /* special case:  d/dA <a|U[A]|b> = - <a|U[A]|db/dB> */
		  (C==B && 0!=shiftb))     /* special case:  d/dB <a|U[B]|b> = - <da/dA|U[B]|b> */
		continue;

	      /* FIXME - type1 integral */
	      /* if ((shell_radius[A][s_A] + shell_radius[B][s_B]) < d_AB) */
	      /* 	continue; */

	      if (t1->smallGrid->start < t1->smallGrid->end) {
		chi = calcChi(t1, U_L, U[C], rAC, dAC, rBC, dBC,
			      la+shifta, sIdx1, pIdx1, shifta, uspA,
			      lb+shiftb, sIdx2, pIdx2, shiftb, uspB);
		if (NULL==chi.array) {
		  result = 1;
		  goto EndOfFunction;
		}
		cb(A, s1, la, shifta, B, s2, lb, shiftb, C, chi.array, args);
		doubleArray_free(chi);
	      }	      
	      
	      /* type2 integral */
	      if (t2->smallGrid->start < t2->smallGrid->end) {
		gamma = calcGamma(t2, FTab, UTab, U[C], dAC, dBC, 
				  la+shifta, sIdx1, pIdx1, omegaA, uspA, shifta, 
				  lb+shiftb, sIdx2, pIdx2, omegaB, uspB, shiftb);		
		if (NULL==gamma.array) {
		  result = 2;
		  goto EndOfFunction;
		}
		cb(A, s1, la, shifta, B, s2, lb, shiftb, C, gamma.array, args);
		doubleArray_free(gamma);
	      }
	    }
	    
	    pIdx2 += h->contraction[sIdx2];
	    sIdx2++;
	  }
	  pIdx1 += h->contraction[sIdx1];
	  sIdx1++;
	}
	sOffset2 = sIdx2;
	pOffset2 = pIdx2;
	doubleArray_free(omegaB);
	doubleArray_free(uspB);
      }
      sOffset1 = sIdx1;
      pOffset1 = pIdx1;
      doubleArray_free(omegaA);
      doubleArray_free(uspA);
    }
    doubleArray_free(UTab);
    doubleArray_free(FTab);
    free(sg);
  }

 EndOfFunction:
  doubleArray_free(omega);
  free(U_L);

  return result;
}


void
libECP_free(libECPHandle *h) {
  int i;
  if (NULL!=h->U) {
    for(i=0; i<h->nrAtoms; i++) {
      if (NULL!=h->U[i])
	ECP_free(h->U[i]);
    }
    free(h->U);
  }
  if (NULL!=h->maxLAtom)
    free(h->maxLAtom);
  if (NULL!=h->t1) 
    Type1_free(h->t1);
  if (NULL!=h->t2) 
    Type2_free(h->t2);
  if (NULL!=h->bessel) {
    freeBesselFunction(h->bessel);
    free(h->bessel);
  }
  if (NULL!=h->largeGrid)
    GCIntegrationTable_free(h->largeGrid);
  if (NULL!=h->fac)
    free(h->fac);
  if (NULL!=h->dfac)
    free(h->dfac);
  if (0!=h->freeIJK) {
    if(NULL!=h->ijk)
      free(h->ijk);
  }
  if(NULL!=h->ijkIndex)
    free(h->ijkIndex);
  if(NULL!=h->cart2sph)
    free(h->cart2sph);
  if(NULL!=h->poly2sph)
    free(h->poly2sph);
  free(h);
}
