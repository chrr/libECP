/* Copyright (c) 2012, Christoph Reimann */

#include <stdlib.h>
#include <math.h>

#include "bessel.h"

/* calculate bessel function K(z)
   K_l(z) = exp(-z) M_l(z), M modified spherical Bessel function of the first kind
   
   function values are obtained corresponding to 
   R. Flores-Moreno et al., J. Comput. Chem. 27 (2005), 1009. (Appendix B)

   parameters:  lMax - should be (maxL + maxAlpha + 6)
                N - number of abscissas (recommendation: 16*100)
		cutoff - cutoff for series expansion (recommendation: 200)
                    
   possible error: series expansion does not converge */
int
tabBesselFunction(BesselFunctionTable *bessel, int lMax, int N, int cutoff, double accuracy) {
  int i, j, l, m;
  double z, f, s;
  /* power series expansion for bessel function K:
     K_l(z) \approx z^l sum_j F_l(z) / G_l(z)
     F_l(z) = e^(-z) (z^2/2)^j / j
     G_l(z) = (2j+1)!!
     K, F and G will be calculated recursively */
  double F[cutoff+1];
  double G[cutoff+lMax+2];
  int dim = N+1;
  double *K = calloc((lMax+1)*dim, sizeof(double)); 
  double *C = calloc(lMax+1, sizeof(double));
  bessel->lMax = lMax;
  bessel->N = N;
  bessel->cutoff = cutoff;
  bessel->K = K;
  bessel->C = C;
  bessel->dimK = dim;

  K[0] = 1.0;
  for(i=1; i<=N; i++) {
    z = i/(N/16.0);
    j = 0;
    f = z*z/2.0;
    F[j] = exp(-z);
    G[j] = 1.0;
    s = F[j]/G[j];
    l = (int)(0.25*sqrt(1.0+16.0*f));
    
    /* series expansion */
    while(s > accuracy || j <= l) {
      K[i] += s;
      j++;
      if (j > cutoff) {
	/* Error: series not converged in bessel function tabulation */
	return 1;
      }
      F[j] = F[j-1]*f/j;
      G[j] = G[j-1]*(2.0*j+1.0);
      s = F[j]/G[j];
    }
    for(l=1; l<=lMax; l++) {
      G[j+l] = G[j+l-1]*(2*j+2*l+1);
    }	      
    
    /* calculate K_l with K_0 expansion */
    f = z;
    for(l=1; l<=lMax; l++) {
      s = 0;
      for(m=0; m<j; m++) 
	s += F[m]/G[l+m];
      K[l*dim+i] = f * s;
      f *= z;
    }
  }

  /* recurrence relation factors */
  for(i=1; i<=lMax; i++) 
    C[i] = i/(2.0*i+1.0);

  return 0;
}


void
freeBesselFunction(BesselFunctionTable *bessel) {
  free(bessel->K);
  free(bessel->C);
}


/* calculate modified spherical Bessel function K(z) weighted with
   an exponential factor
     K(z) := M(z) e^(-z)
   The exponential function e^(-z) is an envelope of the Bessel function
   M(z) for z>=0, so that the K(z) is restricted to the interval [0,1]. 
   
   lit. L.E. McMurchie, E. Davidson, J. Comp. Phys. 44, 289 (1981)
        G. Arfken, H.J. Weber, Mathematical Methods for Physicists,
               Fourth Edition, Academic Press London, 1995, pp. 688 */
void
weightedBesselFunction(const BesselFunctionTable *bessel, const int dim, 
		       const int lmax, const double z, double *K) {
  int i, j, l;
  const double small = 1.0E-7;
  double scale, dz;
    int maxLambda;
  int N = bessel->N;
  int index;
  /* derivatives */
  int ddim = lmax + 6;
  double *dK_i = calloc(ddim, sizeof(double));
  double *dK_j = calloc(ddim, sizeof(double));
  double *A, f;
  ddim = bessel->dimK;
	     
  /* a) zeroth order, z=0 */
  if (z < small) {
    if (z<=0) {
      K[0] = 1.0;
      for(l=1; l<=lmax; l++) 
	K[l*dim] = 0.0;     
    }
    else {
      /* K_l(z) \approx (1-z)z^l/(2l+1)!!
	 => K_0 = 1-z
	    K_1 = (1-z)z/(1*3)
	    K_2 = (1-z)z^2/(1*3*5)
	    K_3 = (1-z)z^3/(1*3*5*7) */
      K[0] = 1-z;
      for(l=1; l<=lmax; l++) 
	K[l*dim] = K[(l-1)*dim]*z/(2*l+1);
    }
  }
  /* b) 0 < z < 16; Taylor series expansion around tabulated function values */
  else if (z < 16.0) {
    maxLambda = lmax + 5;
    scale = N/16;

    /* index of abscissa z in array of tabulated bessel function values */
    index = (int)floor(z*scale+0.5);
    /* dz = z-z_0 */
    dz = z - index/scale;
    scale = 1.0;

    /* initialize dK_i and K */
    for(l=0; l<=lmax; l++) {
      K[l*dim] = dK_i[l] = bessel->K[l*ddim+index];
    }
    for(l=lmax+1; l<=maxLambda; l++) 
      dK_i[l] = bessel->K[l*ddim+index];
    
    /* only 5 terms are necessary to get the taylor expansion of K
       according to R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009--1019. */
    for(i=1; i<=5; i++) {
      index = maxLambda-i;
      
      /* (i-1)th derivatives */
      for(j=0; j<=index+1; j++)
	dK_j[j] = dK_i[j];

      /* ith derivatives */
      dK_i[0] = dK_j[1]-dK_j[0];
      for(j=1; j<=index; j++) {
	/* K_j^(n+1) = j/(2j+1) (K_(j-1)^n+K_(j+1)^n) + 1/(2j+1) K_(j+1)^n */
	/* dK_i[j] = (j/(2.0*j+1.0))*(dK_j[j-1]+dK_j[j+1]) + 1.0/(2.0*j+1.0)*dK_j[j+1] */
	dK_i[j] = bessel->C[j] * (dK_j[j-1] - dK_j[j+1]) - dK_j[j] + dK_j[j+1];
      }
      scale = scale*dz/i;
      for(j=0; j<=lmax; j++)
	K[j*dim] += scale*dK_i[j];
    }
  }
  /* c) z>16: (very accurate) asymptotic formula */
  else {
    /* K_l(z) \approx R_l(-z)/(2z) = \sum_{i=0}^l (l+i)!/(i! (l-i)! (-1)^i (2z)^(i+1))
       => K_0(z) = 0.5/z
          K_1(z) = 0.5/z -  2/(2z)^2
	  K_2(z) = 0.5/z -  6/(2z)^2 + 12/(2z)^3
	  K_3(z) = 0.5/z - 12/(2z)^2 + 60/(2z)^3 - 120/(2z)^4 */
    A = calloc(lmax+1, sizeof(double));
    A[0] = 0.5/z;
    for(l=0; l<=lmax; l++)
      K[l*dim] = A[0];
    for(l=1; l<=lmax; l++) {
      f = l*(l+1);
      for(i=1; i<l; i++) {
	K[l*dim] += f*A[i];
	f *= (l+i+1)*(l-i);
      }
      A[l] = -A[0]*A[l-1]/l;
      K[l*dim] += f*A[l];
    }
    free(A);
  }
  
  free(dK_i); 
  free(dK_j);
}

