
#include "transformations.h"
#include "dimensions.h"
#include "factorial.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* transformation matrices in this file are - in principle - 5 dimensional,
   but in practice reduce to 2 dimensions:
  
   T [l] [m] [i] [j] [k]
     '--.--' '----.----'
      L_DIM     C_DIM

   spherical dimensions: LM_INDEX(l,m), l=0..lmax
   cartesian dimensions: CIJK_INDEX(l,c), c=0..IJK_DIM(l) */

/* generate a matrix for the transformation of cartesian to spherical orbitals
   ref.: H. Schlegel, M. Frisch, Int. J. Quant. Chem. 54 (1995), 83 (formula 15)

   matrix is contructed from submatrices of the form LxC (one submatrix for each l value)
   arguments: lmax - maximum angular momentum
              xyz - exponents for cartesian GTOs
	      fac - precomputed factorials in the range [0..2*lmax]
   returns: dynamically allocated transformation matrix */
double *
TM_cart2sph(const int lmax, const int* const xyz, const double *fac) {
  int l, m, mm;
  int lx, ly, lz;
  int ijk, c, i, j, k, exponent;
  double s, s1, s2;
  double *TM, *T;
  int cnt = 0;
  T = TM = calloc(TM_DIM(lmax,fac), sizeof(double));
  
  for(l=0; l<=lmax; l++) {
    for(m=-l; m<=l; m++) {
      mm = abs(m);
      
      for(c=0; c<IJK_DIM(l); c++) {
	ijk = CIJK_INDEX(l,c);
	/* l = lx+ly+lz */
	lx = xyz[ijk];
	ly = xyz[ijk+1];
	lz = xyz[ijk+2];      
	  
	j = lx+ly-mm;
	if (j<0 || j%2==1) *T = 0.0;
	else {
	  j = j/2;
	  s1 = 0.0;

	  for(i=0; i<=(l-mm)/2; i++) {
	    s2 = 0.0;
	    for(k=0; k<=j; k++) {
	      s = 0.0;
	      if( (m<0 && abs(mm-lx)%2==1) || 
		  (m>0 && abs(mm-lx)%2==0) ) {
		exponent = (mm-lx+2*k)/2;
		/* factor sqrt(2.0) arises because matrix elements are
		   calculated for real spherical harmonics rather than
		   complex ones */
		s = pow(-1.0, exponent)*sqrt(2.0);
	      }
	      else if (m==0 && lx%2==0) {
		exponent = -lx/2+k;
		s = pow(-1.0, exponent);
	      }
	      s2 += n_over_k(j, k, fac)*
		n_over_k(mm,(lx-2*k), fac)*s;
	    }
	    s1 += n_over_k(l,i,fac)*n_over_k(i,j,fac)*
	      pow(-1.0,i)*fac[2*l-2*i]/fac[l-mm-2*i]*s2;
	  }
	  *T = sqrt((fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l]*fac[l-mm]) / 
		     (fac[2*l]*fac[lx]*fac[ly]*fac[lz]*fac[l+mm]))*
                1 / (pow(2.0, l)*fac[l]) * s1;
	}
	T++; cnt++;
      }
    }
  }

  return TM;
}

/* generate a matrix for the transformation from spherical to  cartesian orbitals
   ref.: H. Schlegel, M. Frisch, Int. J. Comput. Chem. 54 (1995), 83 (formula 15) */
double *
TM_sph2cart(double *M, int lmax, const int* const xyz, const double *fac) {
  int l, m;
  int c1, c2, ijk;
  int lx, ly, lz;
  int lx1, ly1, lz1;
  int lx2, ly2, lz2;
  double s, s1, s2;
  int dim = 0;
  double *TM, *T;
  for(l=0; l<=lmax; l++) {
    dim += TM_DIM(l,fac);
  }
  T = TM = calloc(dim, sizeof(double));

  for(l=0; l<=lmax; l++) {
    for(c1=0; c1<IJK_DIM(l); c1++) {
      ijk = CIJK_INDEX(l,c1);
      lx1 = xyz[ijk];
      ly1 = xyz[ijk+1];
      lz1 = xyz[ijk+2];
      s1 = sqrt((fac[lx1] * fac[ly1] * fac[lz1])/
		(fac[2*lx1] * fac[2*ly1] * fac[2*lz1]));

      for(c2=0; c2<IJK_DIM(l); c2++) {
	ijk = CIJK_INDEX(l,c2);
	lx2 = xyz[ijk];
	ly2 = xyz[ijk+1];
	lz2 = xyz[ijk+2];

	/* cartesian components */
	lx = lx1 + lx2;
	ly = ly1 + ly2;
	lz = lz1 + lz2;

	if ((lx%2==0) && (ly%2==0) && (lz%2==0)) {
	  s2 = sqrt((fac[lx2] * fac[ly2] * fac[lz2]) /
		    (fac[2*lx2] * fac[2*ly2] * fac[2*lz2]));
	  s = fac[lx] * fac[ly] * fac[lz] * s1 * s2 /
              (fac[lx/2] * fac[ly/2] * fac[lz/2]);
	  
	  for(m=-l; m<=l; m++) 
	    *T += s * (*M);	  
	}
	T++;
	M++;
      }
    }
  }
  return TM;
}


/* coefficients for the expansion of (unitary sphere) polynomials in terms of real spherical harmonics
   acc. to R. Flores-Moreno et al. J. Comput. Chem. 27 (2006), 1009--1019. (eq. 36) */
double *
TM_poly2sph(double *M, int lmax, const int* const xyz, const double *fac2) {
  unsigned int c1, c2, ijk;
  unsigned int l1, l2, lx1, ly1, lz1, m;
  unsigned int lx2, ly2, lz2;
  unsigned int lx, ly, lz;
  double s, s1, s2;
  double sum;
  unsigned int cdim = C_DIM(lmax);
  unsigned int ldim = L_DIM(lmax);
  double *TM, *T, *TCS;
  
  T = TM = calloc(cdim*ldim, sizeof(double));
    
  for(l1=0; l1<=lmax; l1++) {
    for(c1=0; c1<IJK_DIM(l1); c1++) {
      ijk = CIJK_INDEX(l1,c1);
      lx1 = xyz[ijk];
      ly1 = xyz[ijk+1];
      lz1 = xyz[ijk+2];
      T = TM + (C_DIM(l1-1)+c1)*ldim;
      TCS = M;

      for(l2=0; l2<=l1; l2++) {
	s1 = 4.0 * M_PI * fac2[2*l2+1];
	for(m=0; m<2*l2+1; m++) {
	  sum = 0.0;
	  cdim = IJK_DIM(l2);

	  for(c2=0; c2<cdim; c2++) {
	    ijk = CIJK_INDEX(l2,c2);
	    lx2 = xyz[ijk];
	    ly2 = xyz[ijk+1];
	    lz2 = xyz[ijk+2];
	    
	    lx = lx1+lx2;
	    ly = ly1+ly2;
	    lz = lz1+lz2;
	    
	    if( (lx%2 == 0) &&
		(ly%2 == 0) &&
		(lz%2 == 0) ) {
	      s = s1; 
	      s2 = 1.0/fac2[l1+l2+1];
	      if (lx>2) s2 *= fac2[lx-1];
	      if (ly>2) s2 *= fac2[ly-1]; 
	      if (lz>2) s2 *= fac2[lz-1]; 
	      if (lx2>1) s /= fac2[2*lx2-1]; 
	      if (ly2>1) s /= fac2[2*ly2-1]; 
	      if (lz2>1) s /= fac2[2*lz2-1];
	      sum += sqrt(s) * s2 * (*TCS);
	    }
	    TCS++;
	  }
	  (*T) = sum;
	  T++;
	}
      }
    }
  }
  return TM;
}
