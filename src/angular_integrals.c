/* Copyright (c) 2012, Christoph Reimann */

#include <stdlib.h>
#include <math.h>

#include "angular_integrals.h"
#include "transformations.h"
#include "dimensions.h"
#include "spherical_harmonics.h"


/* tabulate the angular integral Omega(l,m,lambda,mu,i,j,k)
   which is independent of any particular basis set and ECP parameter
   all integrals within the range [0,(L+lmax(basis)] are calculated */
doubleArray
tabOmega(Type2 *t2) {
  int l, m, lambda, mu;
  int alpha, minAlpha;
  int ax, ay, az, d;
  int i, j, k, c, p;
  double *fac = t2->fac;
  double *dfac = t2->dfac;
  double *cijk;
  /* special transformation matrices (initialized during type2 startup) */
  double *T = t2->cart2sph;
  double *P = t2->poly2sph;
  const int tmDim = L_DIM(t2->tmDim);
  /* array increments */
  doubleArray Omega = doubleArray_new(3, L_DIM(t2->maxL), L_DIM(t2->maxLambda), C_DIM(t2->maxAlpha));
  const int inc1 = Omega.dim[3]; 
  const int inc2 = Omega.dim[2]*Omega.dim[3]; 
  /* normalization constant */
  double N;
  const int *ijk = t2->ijk;
  const int *ijkIndex = t2->ijkIndex;
  const int ijkDim = t2->ijkDim;

  for(lambda=0; lambda<=t2->maxLambda; lambda++) {
    for(l=0; l<t2->maxL; l++) {
      i = (lambda+l)%2; 
      j = lambda-l;
      minAlpha = (i>j) ? i : j;

      for(alpha=minAlpha; alpha<=t2->maxAlpha; alpha+=2) {
	for(mu=0; mu<2*lambda+1; mu++) {
	  for(m=0; m<2*l+1; m++) {
	    /* get submatrix for current l,m and lambda,mu pairs */
	    cijk = Omega.array + LM_INDEX(l,m)*inc2 + LM_INDEX(lambda,mu)*inc1 + C_DIM(alpha-1);
	    
	    /* evaluate cartesian components */
	    for(c=0; c<IJK_DIM(alpha); c++) {
	      p = CIJK_INDEX(alpha,c);
	      i = ijk[p];
	      j = ijk[p+1];
	      k = ijk[p+2];

	      /* alpha == 0 => use orthonormality */
	      if(0 == alpha) {
		if(l==lambda && m==mu)
		  cijk[c] = 1.0;
	      }
	      else {
		if(lambda <= l+alpha) {
		  for(d=0; d<IJK_DIM(l); d++) {
		    p = CIJK_INDEX(l,d);
		    ax=ijk[p];
		    ay=ijk[p+1];
		    az=ijk[p+2];

		    /* new normalization factor */
		    N = 0.25*dfac[2*l+1]/M_PI;
		    if (ax>1)
		      N /= dfac[2*ax-1];
		    if (ay>1)
		      N /= dfac[2*ay-1];
		    if (az>1)
		      N /= dfac[2*az-1];
		    N = sqrt(N);
		    		    
		    ax += i;
		    ay += j;
		    az += k;
		    /* determine index for P transformation matrix
		       first calculate index for ijk_index array 
		       (needed to get the index corresponding to (ax+i, ay+j, az+k)) */
		    p = ax*ijkDim*ijkDim + ay*ijkDim + az;
		    /* final index for P matrix */
		    p = ijkIndex[p]*tmDim + LM_INDEX(lambda,mu);
		    
		    cijk[c] += N * T[TM_INDEX(l,m,d,fac)] * P[p];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return Omega;
}


doubleArray 
evalAngularIntegrals(Type2 *t2, doubleArray Omega, int lmax_a, double *rA) {
  int lambda, mu, l, m, c;
  double R[3];
  double *rsph;
  double ijk;
  int la, lmax = t2->maxL - 1 + lmax_a;
  /* array increments for Omega matrix (tabulated angular integral values) */
  const int inc1 = Omega.dim[3]; 
  const int inc2 = Omega.dim[2]*Omega.dim[3];
  /* array increments for A matrix ( A(lambda,l,m,i,j,k) = sum_mu S(lambda,mu) Omega(l,m,lambda,mu,i,j,k) ) */
  doubleArray ang = doubleArray_new(3, lmax+1, L_DIM(t2->maxLambda), C_DIM(lmax_a));
  const int incA1 = ang.dim[3];
  const int incA2 = ang.dim[2]*ang.dim[3];
  double *A = ang.array;

  sphericalCoordinates(rA, R);
  rsph = realSphericalHarmonics(lmax, R[1], R[2], t2->fac, t2->dfac);

  for(l=0; l<t2->maxL; l++) {
    for(m=0; m<2*l+1; m++) {
      for(lambda=0; lambda<=lmax; lambda++) {
	for(la=0; la<=lmax_a; la++) {
	  for(c=0; c<IJK_DIM(la); c++) {
	    ijk = 0.0;
	    for(mu=0; mu<2*lambda+1; mu++) {
	      ijk += rsph[LM_INDEX(lambda,mu)]*
		Omega.array[LM_INDEX(l,m)*inc2 + LM_INDEX(lambda,mu)*inc1 + C_INDEX(la,c)];
	    }
	    A[lambda*incA2 + LM_INDEX(l,m)*incA1 + C_INDEX(la,c)] = ijk;
	  }
	}
      }
    }
  }

  free(rsph);  
  return ang;
}

