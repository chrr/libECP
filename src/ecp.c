
#include <stdlib.h>
#include <math.h>

#include "ecp.h"


int
checkECP(const ECP *U) {
  /* make sure that U->U_l is sorted according to ang. mom. */
  int i, l=0;
  ECPGauss *U_l = U->U_l;
  for (i=0; i<U->N; i++) {
    if (l > U_l->l || U->L < U_l->l) {
      /* not well-ordered or invalid ang. mom. */
      return 1;
    }
    else if (l<U_l->l)
      l = U_l->l;
    U_l++;
  }
  return 0;
}


ECP *
newECP(const int L, const int N) {
  int i;
  ECP *U = malloc(sizeof(ECP));
  U->L = L;
  U->N = N;
  U->U_l = calloc(N, sizeof(ECPGauss));
  for (i=0; i<N; i++)
    U->U_l[i].l = -1;
  return U;
}


/* calculate the radial potential U_l at a grid point r */
double 
evalECP(const ECP *U, const int l, const double r) {
  int i;
  double a, d, n, potential = 0.0;
  const double r2 = r*r;
  ECPGauss *U_l = U->U_l;
  
  for(i=0; i<U->N; i++) {    
    /* loop over shells with angular momentum l */
    if (l == U_l->l) {
      a = U_l->a;
      d = U_l->d;
      n = U_l->n;
      potential += pow(r, n) * d * exp(-a*r2);
    }
    U_l++;
  }
  
  return potential;
}


void
ECP_free(ECP *E) {
  free(E->U_l);
  free(E);
}
