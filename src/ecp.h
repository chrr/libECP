
#ifndef ECP_H
#define ECP_H 1

typedef struct _ECPGauss ECPGauss;
typedef struct _ECP ECP;

struct _ECPGauss {
  int l;
  double a, d, n;
};

struct _ECP {
  /* pseudopotential expanded in Gaussians */
  ECPGauss *U_l;
  /* no. of Gaussians */
  int N;
  /* max. angular momentum */
  int L;
};

/* allocate memory dynamically for new ECP structure */
ECP * newECP(const int L, const int N);

/* check whether ECP structure is well-ordered and angular momentum is in range [0,L] */
int checkECP(const ECP *U);

/* get Gaussian expansion of l-dependent ECP operator U_l 
   parameters: U - ECP
               l - angular momentum 
	       U_l - automatically allocated array of Gaussians representing U_l
   returns: number of Gaussians in U_l */
int getECP(const ECP *U, const int l, ECPGauss **U_l);

/* evaluate the l-dependent ECP operator U_l at point r */
double evalECP(const ECP *U, const int l, const double r);

void ECP_free(ECP *E);

#endif
