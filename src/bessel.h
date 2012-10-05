
#ifndef _BESSEL_H
#define _BESSEL_H 1

/* modified spherical Bessel function of the first kind */
typedef struct {
  /* lMax: max. value for lambda */
  int lMax;
  /* N: number of abscissas */
  int N;
  /* series expansion cutoff */
  int cutoff;
  int lambdaMax;
  /* row dimension of K */
  int dimK;
  double *K;
  /* coefficients in recurrence relation for derivatives of bessel functions */
  double *C;
} BesselFunctionTable;

/* tabulate exponentially weighted bessel function K_l(z) = exp(-z) M_l(z)
   M: modified spherical Bessel function of the first kind

   function values are obtained corresponding to 
   R. Flores-Moreno et al., J. Comput. Chem. 27 (2005), 1009. (Appendix B)

   parameters:  K - array (0..lMax)x(0..N) to hold the bessel function values */
int tabBesselFunction(BesselFunctionTable *bessel, int lMax, int N, int cutoff, double accuracy);

/* calculate exponentially weighted modified spherical Bessel function K(z) 

   lit. L.E. McMurchie, E. Davidson, J. Comp. Phys. 44, 289 (1981)
          G. Arfken, H.J. Weber, Mathematical Methods for Physicists,
          Fourth Edition, Academic Press London, 1995, pp. 688 */
void weightedBesselFunction(const BesselFunctionTable *bessel, const int dim, const int lmax, const double z, double *K);

void freeBesselFunction(BesselFunctionTable *bessel);

#endif /* _BESSEL_H */
