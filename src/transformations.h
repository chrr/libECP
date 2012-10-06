/* Copyright (c) 2012, Christoph Reimann */

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H 1

/* transformation submatrix dimensions (M_DIM x IJK_DIM):
   l:      0    1    2    3    4    5    6
   MxIJK:  1    9   30   70  135  231  364
   TM_DIM: 1   10   40  110  245  476  

   => LxC:    (2l^3 + 7l^2 + 7l + 2)/2 = (l+1/2)(l+1)(l+2)
      TM_DIM: (3*(l+1)-1)*factorial((l+1)+2)/(12*factorial((l+1)-1)), l=0..lmax */
#define TM_DIM(l,fac) ((int)((3*(l)+2)*fac[l+3]/(12*fac[l])))
#define TM_INDEX(l,m,c,fac) ((l)==0 ? 0 : (TM_DIM(l-1,fac)+(m)*IJK_DIM(l)+c))

/* generate a matrix for the transformation of cartesian to spherical orbitals
   ref.: H. Schlegel, M. Frisch, Int. J. Quant. Chem. 54 (1995), 83 (formula 15)
   arguments: lmax - maximum angular momentum
              xyz - exponents for cartesian GTOs
	      fac - precomputed factorials in the range [0..2*lmax]
   returns: dynamically allocated transformation matrix */
double * TM_cart2sph(int lmax, const int *xyz, const double *fac);
double * TM_sph2cart(double *M, int lmax, const int* const xyz, const double *fac);

/* coefficients for the expansion of (unitary sphere) polynomials in terms of real spherical harmonics
   acc. to R. Flores-Moreno et al. J. Comput. Chem. 27 (2006), 1009--1019. (eq. 36)
   returns: dynamically allocated C_DIM(lmax) x L_DIM(lmax) transformation matrix T  
   dimensions:
     l:        0    1    2    3    4
     C x:     1x   4x  10x  20x  35x  L_DIM(lmax)
     TM_DIM = C_DIM(lmax) x L_DIM(lmax) */   
double * TM_poly2sph(double *M, int lmax, const int* const xyz, const double *fac2);

#endif
