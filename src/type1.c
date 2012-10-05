
#include <math.h>

#include "ecp_array.h"
#include "type1.h"
#include "ecp.h"
#include "gc_integrators.h"
#include "dimensions.h"
#include "spherical_harmonics.h"
#include "util.h"


Type1 *
Type1_new(int nrAtoms, double *geometry, ECP **U, 
	  int *shells, int *am, int *contraction, double *d, double *a, 
	  int derivative, const double tol) {
  Type1 *t1 = malloc(sizeof(Type1));
  t1->contraction = contraction;
  t1->d = d;
  t1->a = a;
  t1->shift = derivative;
  t1->smallOrder = 128;
  t1->tolerance = tol;

  return t1;
}


void
Type1_free(Type1 *t1) {
  GCIntegrationTable_free(t1->smallGrid);
  free(t1);
}


void
Type1_init(Type1 *t1, const int maxLambda, double *fac, double *dfac,
	   int *ijk, int *ijkIndex, const int ijkDim,
	   double *poly2sph, const int tmDim, BesselFunctionTable *bessel,
	   int largeOrder, GCIntegrationTable *largeGrid, 
	   const double accuracy, const double lnAccuracy) {
  t1->maxLambda = maxLambda;  
  t1->fac = fac;
  t1->dfac = dfac;
  t1->ijk = ijk;
  t1->ijkIndex = ijkIndex;
  t1->ijkDim = ijkDim;
  t1->poly2sph = poly2sph;
  t1->tmDim = tmDim;
  t1->bessel = bessel;
  t1->accuracy = accuracy;
  t1->lnAccuracy = lnAccuracy;
  /* grids for radial integration */
  t1->smallGrid = initGridGC_PS93(t1->smallOrder, t1->tolerance);  
  /* transform small grid from (-1,+1) to (0,\infty) */
  transformGrid_KK(t1->smallGrid->order, t1->smallGrid->x, t1->smallGrid->w);
  t1->largeOrder = largeOrder;
  t1->largeGrid = largeGrid; 
}


typedef struct {
  double C;
  double *exponents;
  /* weighted Bessel functions K_lambda(z) */
  double *K;
  /* tabulated effective core potential */
  double *U;
  /* r^(alpha+beta) */
  double *rn;
  /* accuracy w.r.t. the exponential function */
  double minExp;
} QIntegrandParameters;


double
QIntegrand(double r, unsigned int index, const QIntegrandParameters *p) {
  double Q = 0.0;
  double expo = p->exponents[index];

  if(expo >= p->minExp) {
    Q = p->C * p->rn[index] * p->U[index] * p->K[index] * exp(expo);
  }

  return Q;
}


/* alternative: calculate \sum_l1 \sum_l2 T(l1,l2,lab) for a combination of shells with the ranges
     l1,l2: 0..L
   by integration over primitive Gaussian functions */
doubleArray
calcQ(Type1 *t1, double *U_L, ECP *U, const double s, const int lab, 
      const double dAC, double ca, const double za, const double shifta,
      const double dBC, double cb, const double zb, const double shiftb) {
  GCIntegrationTable *grid;
  QIntegrandParameters ps_Q;
  Integrand f_Q = { (IntegrandFunction)QIntegrand, &ps_Q };
  int l1, l2, n;
  double zd2 = -za*dAC*dAC-zb*dBC*dBC;
  double zeta_p, p, r;
  double *K, *rn;
  doubleArray integrals = doubleArray_new(2, lab+1, lab+1);
  double *T = integrals.array;
  const int Tinc = integrals.dim[2];
  double z = -za-zb;
  int nFailed = 0;
  int *indicesFailed = calloc(2*(lab+1)*(lab+1), sizeof(int));

  ps_Q.C = ca*cb*exp(zd2);
  ps_Q.exponents = calloc(t1->largeGrid->order, sizeof(double));
  ps_Q.minExp = t1->lnAccuracy;
  ps_Q.U = calloc(t1->largeGrid->order, sizeof(double));
  K = calloc((lab+1)*t1->largeGrid->order, sizeof(double));
  rn = calloc((lab+1)*t1->largeGrid->order, sizeof(double));
  
  /* 1. try to evaluate integral using the small generic grid */
  for (n=t1->smallGrid->start; n<t1->smallGrid->end; n++) {
    r = t1->smallGrid->x[n];
    ps_Q.exponents[n] = (z*r + s)*r;
    ps_Q.U[n] = U_L[n];
    weightedBesselFunction(t1->bessel, t1->smallGrid->order, lab, s*r, K + n);
    rn[n] = 1.0;
    for(l1=1; l1<=lab; l1++) {
      rn[l1*t1->smallGrid->order + n] = r * rn[(l1-1)*t1->smallGrid->order + n];
    } 
  }
  /* perform numerical integration with PS93 two-point sequence */      
  for(l1=0; l1<=lab; l1++) {
    ps_Q.rn = rn + l1*t1->smallGrid->order;
    /* angular integral is non-zero only if alpha+beta-lambda is even */
    for(l2=l1; l2>=0; l2-=2) {
      ps_Q.K = K + l2*t1->smallGrid->order;
      if (0 != integrateGC_PS93(&f_Q, t1->smallGrid)) {
	indicesFailed[nFailed*2] = l1;
	indicesFailed[nFailed*2+1] = l2;
	nFailed++;
      }
      else {
	T[l1*Tinc + l2] += t1->smallGrid->I;
      }
    }
  }

  /* 2. in case of failure, use large grid with linear grid point mapping */
  if (0<nFailed) {
    /* regroup exponential terms together due to numerical issues */
    ps_Q.C = ca*cb;
    grid = GCIntegrationTable_copy(t1->largeGrid);
    zeta_p = za + zb;
    p = (za*dAC + zb*dBC)/zeta_p;
    transformGrid_FM06(grid->order, grid->x, grid->w, zeta_p, p);
    
    /* tabulate r^(alpha+beta), bessel functions K(r), exp(), U_L(r) */
    for(n=0; n<grid->order; n++) {
      r = grid->x[n];
      /* formula for Q exponent modified as we are using weighted bessel functions 
	 cmp. Flores-Moreno et al., JCP 27 (2006), 1009--1019, eq. (24) */
      ps_Q.exponents[n] = (z*r + s)*r + zd2;
      
      /* simple screening: check for exp(t1->lnAccuracy) > 0.0 */
      if(r>dAC && r>dBC && ps_Q.exponents[n]<ps_Q.minExp) {
	grid->end = n-1;
	break;
      }
      else if(ps_Q.exponents[n]>=ps_Q.minExp) {
	weightedBesselFunction(t1->bessel, grid->order, lab, s*r, K + n);
	ps_Q.U[n] = evalECP(U, U->L, r);
	/* FIXME - when calculating r^(alpha+beta):
	   what is the convention for the exponent here - alpha+beta+2 or alpha+beta?
	   depends on whether the +2 has been included in the parameter set or not... */
	rn[n] = 1.0;
	for(l1=1; l1<=lab; l1++) {
	  rn[l1*grid->order + n] = r * rn[(l1-1)*grid->order + n];
	}
      }
    }

    /* perform numerical integration with PSM92 sequence */
    for(n=0; n<nFailed; n++) {
      l1 = indicesFailed[n*2];
      l2 = indicesFailed[n*2+1];
      ps_Q.rn = rn + l1*grid->order;
      ps_Q.K = K + l2*grid->order;
      if (0 != integrateGC_PSM92(&f_Q, grid)) {
	doubleArray_free(integrals);
	goto endOfFunction;
      }
      else {
	T[l1*Tinc + l2] += grid->I;
      }
    }
  }
  
 endOfFunction:
  free(ps_Q.exponents);
  free(ps_Q.U);
  free(K);
  free(rn);
  free(indicesFailed);
  if (0<nFailed) {
    GCIntegrationTable_free(grid);
  }
  return integrals;
}


doubleArray
calcChi(Type1 *t1, double *U_L, ECP *U,
	double *rAC, const double dAC, double *rBC, const double dBC, 
	const int la, const int shella, const int offseta, const int shifta, doubleArray uspA,  
	const int lb, const int shellb, const int offsetb, const int shiftb, doubleArray uspB) {
  const int N_a = t1->contraction[shella];
  const double *zs_a = t1->a + offseta;
  const double *cs_a = t1->d + offseta;
  const int N_b = t1->contraction[shellb];
  const double *zs_b = t1->a + offsetb;
  const double *cs_b = t1->d + offsetb;
  int pa, pb;
  double ca, za, cb, zb, factor;
  double P[3], S[3];
  const int lab = la+lb;
  int ax, ay, az, bx, by, bz;
  int i, j, p, l, m, lmax, lx, ly, lz;
  double *rsph, *PM = t1->poly2sph;
  int *ijkIndex = t1->ijkIndex;
  const int ijkDim = t1->ijkDim;
  const int tmDim = L_DIM(t1->tmDim);
  doubleArray Q, chi = doubleArray_new(2, C_DIM(la), C_DIM(lb));
  int incQ, inc = chi.dim[2];
  
  for(pa=0; pa<N_a; pa++) {
    za = zs_a[pa];
    ca = cs_a[pa]; 
    /* nth derivative: extra factor of zeta^n */
    for(p=0; p<shifta; p++)
      ca *= za;
    
    for(pb=0; pb<N_b; pb++) {
      zb = zs_b[pb];
      cb = cs_b[pb];
      for(p=0; p<shiftb; p++)
	cb *= zb;

      /* note: use of rAC,rBC introduces change in sign */
      for (p=0; p<3; p++) 
	P[p] = 2.0*(za*rAC[p] + zb*rBC[p]);
      sphericalCoordinates(P, S);
      rsph = realSphericalHarmonics(lab, S[1], S[2], t1->fac, t1->dfac);
            
      /* bessel weight 's' corresponds to S[0] */
      Q = calcQ(t1, U_L, U, S[0], lab, dAC, ca, za, shifta, dBC, cb, zb, shiftb);
      if (NULL==Q.array) {
	doubleArray_free(chi);
	free(rsph);
	return chi;
      }
      incQ = Q.dim[2];

      /* construct cartesian components */
      // FIXME - ??? for(alpha=0; alpha<=la; alpha++) {
      //               for(c1=0; c1<IJK_DIM(alpha); c1++) { ...
      for(ax=0; ax<=la; ax++) {
	for(ay=0; ay<=la-ax; ay++) {
	  for(az=0; az<=la-ax-ay; az++) {
	    i = ijkIndex[ax*ijkDim*ijkDim + ay*ijkDim + az];
	    for(bx=0; bx<=lb; bx++) {
	      lx = ax + bx;
	      for(by=0; by<=lb-bx; by++) {
		ly = ay + by;
		for(bz=0; bz<=lb-bx-by; bz++) {
		  j = ijkIndex[bx*ijkDim*ijkDim + by*ijkDim + bz];
		  lz = az + bz;
		  lmax = lx + ly + lz;

 		  /* final index for P matrix */
		  p = ijkIndex[lx*ijkDim*ijkDim + ly*ijkDim + lz];

		  for(l=lmax; l>=0; l-=2) {
		    factor = 0.0;
		    for(m=0; m<2*l+1; m++) {
		      factor += rsph[LM_INDEX(l,m)] * PM[p*tmDim+LM_INDEX(l,m)];
		    }
		    chi.array[i*inc+j] += factor * Q.array[lmax*incQ+l];
		  }
		}
	      }
	    }
	  }
	}
      }
      doubleArray_free(Q);
      free(rsph);
    }
  }

  return calcPolynomials(chi, t1->ijk, t1->ijkIndex, t1->ijkDim,
			 t1->fac, t1->accuracy, 4.0*M_PI, la, uspA, lb, uspB);
}

