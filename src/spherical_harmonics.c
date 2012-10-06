/* Copyright (c) 2012, Christoph Reimann */

#include <math.h>
#include <stdlib.h>

#include "spherical_harmonics.h"
#include "dimensions.h"

#define RSH(l,m) LM_INDEX((l),M_INDEX(l,m))

/* calculate real spherical harmonics
   
   ref. W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery,
        Numerical Recipes, 2nd Edition, Cambridge University Press */
double *
realSphericalHarmonics(const int lmax, const double theta, const double phi, double *fac, double *dfac) {
  int l, m;
  /* x := argument for legendre polynomial */
  double x;
  /* P := associated legendre polynomials */
  double *P; 
  /* s := sin(m \phi), c := cos(m \phi) */
  double *s, *c;
  double norm0, norm;
  /* rsh := real spherical harmonics */
  double *rsh; 

  /* allocate memory */
  s = calloc(lmax+2, sizeof(double)); 
  c = calloc(lmax+2, sizeof(double));
  rsh = calloc(L_DIM(lmax), sizeof(double));
  P = calloc(L_DIM(lmax), sizeof(double));

  /* calculate associated legendre polynomials P(x) */
  x = cos(theta);
  /* special case: theta == 0 */
  if (1.0 == x) {
    for(l=0; l<=lmax; l++) 
      P[RSH(l,0)] = 1.0;
  }
  /* special case: theta == pi */
  else if (-1.0 == x) {
    P[RSH(0,0)] = 1.0;
    for(l=1; l<=lmax; l++) 
      P[RSH(l,0)] = -P[RSH(l-1,0)];
  }
  else {
    s[1] = sqrt(1.0-x*x);
    for(l=2; l<=lmax; l++) 
      s[l] = s[l-1]*s[1];
    for(l=0; l<=lmax; l++) {
      m = l; 
      if (0 == m) {
	P[RSH(l,0)] = 1.0;
      }
      else {
	P[RSH(l,m)] = s[m]*dfac[2*m-1];
	m = l-1;
	P[RSH(l,m)] = x*(2*m+1)*P[RSH(l-1,m)];
	if (l > 1) {
	  for(m=0; m<=l-2; m++) {
	    P[RSH(l,m)] = (x*(2*l-1)*P[RSH(l-1,m)] - 
			   (l+m-1)*P[RSH(l-2,m)]) / (l-m);
	  }
	}
      }
    }
  }

  /* normalization */
  for(l=0; l<=lmax; l++) {
    norm0 = sqrt((2.0*l+1.0)/(2.0*M_PI));
    P[RSH(l,0)] = norm0 * P[RSH(l,0)];
    for(m=1; m<=l; m++) {
      norm = sqrt(fac[l-m]/fac[l+m]) * norm0;
      P[RSH(l,m)] = norm * P[RSH(l,m)];
    }
  }

  /* trigonometric functions */
  if (lmax > 0) {
    if (0.0 == phi) {
      for(m=0; m<=lmax; m++) {
	s[m] = 0.0;
	c[m] = 1.0;
      }
    }
    else {
      s[1] = sin(phi);
      c[1] = cos(phi);
      for(m=2; m<=lmax; m++) {
	s[m] = s[1]*c[m-1] + c[1]*s[m-1];
	c[m] = c[1]*c[m-1] - s[1]*s[m-1];
      }
    }
  }

  /* normalized real spherical harmonics */
  for(l=0; l<=lmax; l++) {
    rsh[RSH(l,0)] = P[RSH(l,0)]/sqrt(2.0);

    for(m=1; m<=l; m++) {
      rsh[RSH(l,-m)] = P[RSH(l,m)]*s[m];
      rsh[RSH(l,+m)] = P[RSH(l,m)]*c[m];
    }
  }

  /* free memory */
  free(P);
  free(s);
  free(c);

  return rsh;
} 
