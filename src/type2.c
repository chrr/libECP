
#include <stdlib.h>
#include <math.h>

#include "type2.h"
#include "util.h"
#include "dimensions.h"
#include "factorial.h"
#include "transformations.h"
#include "angular_integrals.h"

int 
initBS(Type2 *t2, int *shells, int *am, 
       int *contraction, double *d, double *a) {
  int i, lmax = 0;
  t2->shells = shells;
  t2->nrShells = 0;
  t2->maxShells = 0;
  t2->am = am;
  t2->contraction = contraction;
  t2->d = d;
  t2->a = a;
  if (NULL != shells) {
    /* determine total number of shells and max. number of shells per atom */
    for(i=0; i<t2->nrAtoms; i++) {
      if (shells[i] > t2->maxShells)
	t2->maxShells = shells[i];
      t2->nrShells += shells[i];
    }
  }
  /* determine max. angular momentum */
  for(i=0; i<t2->nrShells; i++) {
    if(am[i] > lmax)
      lmax = am[i];
  }  
  return lmax;
}


Type2 *
Type2_new(const int nrAtoms, double *geometry, ECP **U, 
	  int *shells, int *am, 
	  int *contraction, double *d, double *a, 
	  const int derivative, const double tol) {
  int i, maxLECP = 0;
  Type2 *t2 = malloc(sizeof(Type2));
  t2->nrAtoms = nrAtoms;
  t2->geometry = geometry;
  t2->U = *U;
  t2->maxLBS = initBS(t2, shells, am, contraction, d, a);
  for(i=0; i<nrAtoms; i++) {
    if (NULL!=U[i] && U[i]->L>maxLECP)
      maxLECP = U[i]->L;
  }
  t2->maxL = maxLECP;
  t2->shift = derivative;  
  /* assign default values */
  t2->accuracy = 1.0E-14;
  t2->lnAccuracy = log(t2->accuracy)-2;
  t2->smallOrder = 128;
  t2->largeOrder = 1024;
  /* tolerance in quadrature for radial integrals (FM06: 1.0E-8) */
  t2->tolerance = tol;
  /* parameters for bessel function evaluation */
  t2->N = 16*100;
  t2->cutoff = 200; 

  return t2;
}

void
Type2_init(Type2 *t2, double *fac, double *dfac, int *ijk, int *ijkIndex, const int ijkDim,
	   double *cart2sph, double *poly2sph, const int tmDim, GCIntegrationTable *largeGrid, 
	   BesselFunctionTable *bessel) {
  t2->freeMemory = ECP_FALSE;
  t2->maxAlpha = t2->maxLBS + t2->shift;
  t2->maxLambda = t2->maxL - 1 + t2->maxAlpha;
  /* factorials */
  t2->fac = fac;
  t2->dfac = dfac;
  t2->ijk = ijk;
  t2->ijkIndex = ijkIndex;
  t2->ijkDim = ijkDim;
  t2->tmDim = tmDim;
  t2->cart2sph = cart2sph;
  t2->poly2sph = poly2sph;
  /* grids for radial integration */
  t2->smallGrid = initGridGC_PS93(t2->smallOrder, t2->tolerance);  
  /* transform small grid from (-1,+1) to (0,\infty) */
  transformGrid_KK(t2->smallGrid->order, t2->smallGrid->x, t2->smallGrid->w);
  t2->largeGrid = largeGrid;
  t2->bessel = bessel;
}


ErrType2
Type2_initDefault(Type2 *t2) {
  t2->freeMemory = ECP_TRUE;
  t2->maxAlpha = t2->maxLBS + t2->shift;
  t2->maxLambda = t2->maxL - 1 + t2->maxAlpha;
  t2->lmax = t2->maxL + t2->maxAlpha + 6;
  /* factorials (current limit determined by call to TM_DIM) */
  t2->fac = factorials(2*t2->maxLambda+3);
  t2->dfac = double_factorials(2*t2->maxLambda+1);

  /* FIXME - use of maxLambda correct here, or is maxAlpha sufficient?? */
  t2->ijk = cartesianShellOrder(t2->maxLambda);
  t2->ijkIndex = cartesianShellOrderIndex(t2->maxLambda, t2->ijk);
  t2->ijkDim = t2->maxLambda+1;
  t2->tmDim = t2->maxLambda;
  t2->cart2sph = TM_cart2sph(t2->tmDim, t2->ijk, t2->fac);
  t2->poly2sph = TM_poly2sph(t2->cart2sph, t2->tmDim, t2->ijk, t2->dfac);

  /* grids for radial integration */
  t2->smallGrid = initGridGC_PS93(t2->smallOrder, t2->tolerance);  
  t2->largeGrid = initGridGC_PSM92(t2->largeOrder, t2->tolerance);
  /* transform small grid from (-1,+1) to (0,\infty) */
  transformGrid_KK(t2->smallGrid->order, t2->smallGrid->x, t2->smallGrid->w);

  t2->bessel = malloc(sizeof(BesselFunctionTable));
  if (0 != tabBesselFunction(t2->bessel, t2->lmax, t2->N, t2->cutoff, t2->accuracy)) 
    return BESSEL_NOT_CONVERGED;

  return TYPE2_OK;
}

void 
Type2_free(Type2 *t2) {
  if(ECP_TRUE==t2->freeMemory) {
    free(t2->fac);
    free(t2->dfac);
    free(t2->ijk);
    free(t2->ijkIndex);
    free(t2->cart2sph);
    free(t2->poly2sph);
    GCIntegrationTable_free(t2->largeGrid);
    freeBesselFunction(t2->bessel);
    free(t2->bessel);
  }
  GCIntegrationTable_free(t2->smallGrid);
  free(t2);
}

/* screening by basis set cutoff */
ScreenedGrid
basisSetScreening(GCIntegrationTable *grid, const double shellRadius, const double d_IC) {
  int j;
  const double r_min = d_IC - shellRadius;
  const double r_max = d_IC + shellRadius;
  double *r = grid->x;
  ScreenedGrid s;

  /* determine upper bound for grid points */
  for(j=(int)grid->end; j>=(int)grid->start; j--) {
    if (r[j] < r_min)		
      break;
  }
  
  /* index for first grid point that is not screened */
  s.start = j+1;
  for(j=(int)grid->end; j>=(int)grid->start; j--) {
    if (r[j] <= r_max)
      break;
  }

  /* index for last grid not to be screened */
  if (j<0 || r[j]>r_max)
    s.end = -1;
  else
    s.end = j;
  if (s.end >= s.start) 
    s.skipShell = ECP_FALSE;
  else
    s.skipShell = ECP_TRUE;

  return s;
}


/* calculate r^n U_l(r), n=0..lmax, l=0..L-1 and evaluate potential screening */
doubleArray
tabECP_FM06(ECP *U, GCIntegrationTable *grid, const int lmax, ECPBool useScreening, const double tol) {
  int l, n, lab;
  double *P;
  doubleArray U_tab = doubleArray_new(3, U->L, lmax+1, grid->order);
  const int inc1 = U_tab.dim[3];
  const int inc2 = U_tab.dim[2]*U_tab.dim[3];
  double *U_l = calloc(grid->order, sizeof(double));
  double *r = grid->x;

  /* evaluate U_l for each l at every grid point */
  for (l=0; l<U->L; l++) {
    for(n=0; n<grid->order; n++) {
      U_l[n] = evalECP(U, l, r[n]);
    }
    
    /* evaluate potential cutoff */
    if (ECP_TRUE == useScreening) {
      grid->end = potentialScreening(grid, U_l, tol);
    }

    for(n=grid->start; n<=grid->end; n++) {
      /* note: weight factor w[n] is applied automatically during numerical integration */
      P = U_tab.array + l*inc2 + n;
      P[0] = U_l[n];
      
      /* tabulate r^(alpha+beta) prefactor as well */
      for(lab=1; lab<=lmax; lab++) {
	P[lab*inc1] = P[(lab-1)*inc1] * r[n];
      }
    }
  }

  free(U_l);
  return U_tab;
}


doubleArray
calcF_FM06(Type2 *t2, double *rC, ScreenedGrid **screening) {
  int A, i, n, dim = t2->maxShells;
  double dAC, zeta, da, e;
  int sa, pa, *primitives = calloc(t2->nrAtoms*dim, sizeof(int));
  int lmaxA, lmax = t2->maxL - 1;
  /* t2->shift>0 in case derivatives need to be calculated */
  const int nRuns = 1 + t2->shift;
  int run, order;
  GCIntegrationTable *grid = t2->smallGrid;
  ScreenedGrid *sg = (*screening) = calloc(t2->nrShells, sizeof(ScreenedGrid));
  ECPBool *skipAtom = calloc(t2->nrAtoms, sizeof(ECPBool));
  /* dimension for matrix F */
  doubleArray FMat = doubleArray_new(4, nRuns, t2->nrShells, t2->maxLambda+1, t2->smallGrid->order);
  const int inc1 = FMat.dim[4];
  const int inc2 = FMat.dim[3]*FMat.dim[4]; 
  const int inc3 = FMat.dim[2]*FMat.dim[3]*FMat.dim[4];
  double *F = FMat.array;
  double *geom = t2->geometry;
  double *r = t2->smallGrid->x;
  double *K = calloc(t2->maxLambda+1, sizeof(double));  
  int offset, index;

  /* evaluate basis set screening */
  offset = index = 0;
  for(A=0; A<t2->nrAtoms; A++) {
    dAC = distance(rC, &geom[A*3]);
    skipAtom[A] = ECP_TRUE;
    for (sa=0; sa<t2->shells[A]; sa++) {
      e = shellRadius(t2->contraction[index], t2->am[index], 
		      &t2->d[offset], &t2->a[offset], t2->accuracy);
      sg[index] = basisSetScreening(grid, e, dAC);
      if (ECP_FALSE == sg[index].skipShell)
	skipAtom[A] = ECP_FALSE;
      primitives[A*dim+sa] = offset;
      offset += t2->contraction[index];
      index++;
    }
  }
  
  /* tabulate values for F(lambda) */
  for(run=0; run<nRuns; run++) {
    order = run;
    index = -1; 

    for(A=0; A<t2->nrAtoms; A++) {
      if (ECP_TRUE == skipAtom[A]) {
	index += t2->shells[A];
	continue;
      }
      dAC = distance(rC, &geom[A*3]);

      /* tabulate integrals for each shell */
      for (sa=0; sa<t2->shells[A]; sa++) {
	index++;
	if (ECP_TRUE == sg[index].skipShell)
	  continue;
	lmaxA = lmax + t2->am[index] + order;

	/* sum up the values for the primitives */
	for(pa=0; pa<t2->contraction[index]; pa++) {
	  offset = primitives[A*dim+sa] + pa;
	  zeta = t2->a[offset];
	  da = t2->d[offset];

	  /* in case of derivatives: need d^n */
	  for(n=0; n<order; n++)
	    da *= zeta;
	  
	  /* evaluate the integral on the (screened) grid */
	  for(n=sg[index].start; n<sg[index].end; n++) {
	    weightedBesselFunction(t2->bessel, 1, lmaxA, 2.0*zeta*dAC*r[n], K);
	    e = dAC-r[n];
	    e = exp(-zeta*e*e);
	    for(i=0; i<=lmaxA; i++) {
	      F[run*inc3 + index*inc2 + i*inc1 + n] += da * K[i] * e;
	    }
	  }
	}
      }
    }
  }

  free(skipAtom);
  free(primitives);
  free(K);

  return FMat;
}


typedef struct {
  double *Fa, *Fb;
  double *U_tab;
} TIntegrandParameters;


double
TIntegrand_FM06(double x, int n, const TIntegrandParameters *p) {
  return p->Fa[n] * p->Fb[n] * p->U_tab[n];
}									


/* calculate T(l1,l2,lab) for a combination of shells with the ranges
     l1: 0..(l+la)
     l2: 0..(l+lb)
     lab: 0..(la+lb)
   either by using tabulated values or by integration over primitive Gaussian functions */
typedef struct {
  int N;
  int *l1, *l2, *l3;
} IntegrationFailure;


IntegrationFailure
calcT_fast_FM06(doubleArray integrals, Type2 *t2, const int l, const doubleArray U_tab,  
		const int alpha, const int beta, const doubleArray F, double *Fa, double *Fb) {
  int laC = alpha + l;
  int lbC = beta + l;
  int lab = alpha + beta;
  int l1, l2, l3;
  double *T = integrals.array;
  const int Tinc1 = integrals.dim[3];
  const int Tinc2 = integrals.dim[2]*integrals.dim[3];
  TIntegrandParameters ps_T;
  Integrand f_T = { (IntegrandFunction)TIntegrand_FM06, &ps_T };
  /* Fa and Fb have the same dimensions */
  const int Finc1 = F.dim[4];
  const int Uinc1 = U_tab.dim[3];
  const int Uinc2 = U_tab.dim[2]*U_tab.dim[3]; 
  IntegrationFailure failed;
  const int maxFailed = (laC+1)*(lbC+1)*(lab+1);
  failed.N = 0;
  failed.l1 = calloc(maxFailed, sizeof(int));
  failed.l2 = calloc(maxFailed, sizeof(int));
  failed.l3 = calloc(maxFailed, sizeof(int));

  for(l1=0; l1<=laC; l1++) {
    ps_T.Fa = &Fa[l1*Finc1]; 
    for(l2=0; l2<=lbC; l2++) {
      ps_T.Fb = &Fb[l2*Finc1]; 
      for(l3=0; l3<=lab; l3++) {
	ps_T.U_tab = U_tab.array + l*Uinc2 + l3*Uinc1;	
	
	if (0 != integrateGC_PS93(&f_T, t2->smallGrid)) {
	  /* numeric integration on small grid failed */
	  failed.l1[failed.N] = l1;
	  failed.l2[failed.N] = l2;
	  failed.l3[failed.N] = l3;
	  failed.N++;
	}
	else {
	  T[l1*Tinc2 + l2*Tinc1 + l3] = t2->smallGrid->I;
	}
      }
    }
  }

  return failed;
}

typedef struct {
  double C;
  double *exponents;
  /* weighted Bessel functions K_lambda(z) */
  double *K_a, *K_b;
  /* tabulated effective core potential */
  double *U;
  /* r^(alpha+beta) */
  double *rn;
  /* accuracy w.r.t. the exponential function */
  double minExp;
} QIntegrandParameters;


double
QIntegrand_FM06(double r, int index, const QIntegrandParameters *p) {
  double Q = 0.0;
  double expo = p->exponents[index];

  if(expo >= p->minExp) {
    Q = p->C * p->U[index] * p->rn[index]
      * p->K_a[index] * p->K_b[index] 
      * exp(expo);
  }


  return Q;
}

/* alternative: calculate \sum_l1 \sum_l2 T(l1,l2,lab) for a combination of shells with the ranges
     l1: 0..(l+la)
     l2: 0..(l+lb)
     lab: 0..(la+lb)
     by integration over primitive Gaussian functions */
IntegrationFailure
calcT_FM06(doubleArray integrals, IntegrationFailure repeat, Type2 *t2, ECP *U, const int l,
	   const double dAC, const int alpha, const int shella, const int shellOffseta, const int shifta,  
	   const double dBC, const int beta, const int shellb, const int shellOffsetb, const int shiftb) {
  const int Na = t2->contraction[shella];
  const double *zs_a = t2->a + shellOffseta;
  const double *cs_a = t2->d + shellOffseta;
  const int Nb = t2->contraction[shellb];
  const double *zs_b = t2->a + shellOffsetb;
  const double *cs_b = t2->d + shellOffsetb;
  GCIntegrationTable *grid;
  QIntegrandParameters ps_Q;
  Integrand f_Q = { (IntegrandFunction)QIntegrand_FM06, &ps_Q };
  int pa, pb;
  int l1, l2, l3, n;
  double zeta_a, zeta_b, zeta_p, p;
  double s1, s2, d1, d2, r;
  double *K_a, *K_b, *rn;
  const int laC = alpha+l;
  const int lbC = beta+l;
  const int lab = alpha+beta;
  const int Tinc1 = integrals.dim[3];
  const int Tinc2 = integrals.dim[3]*integrals.dim[2];
  double *T = integrals.array;
  const int nFailed = repeat.N-1;

  ps_Q.exponents = calloc(t2->largeGrid->order, sizeof(double));
  ps_Q.minExp = t2->lnAccuracy;
  ps_Q.U = calloc(t2->largeGrid->order, sizeof(double));
  K_a = calloc((laC+1)*t2->largeGrid->order, sizeof(double));
  K_b = calloc((lbC+1)*t2->largeGrid->order, sizeof(double));
  rn = calloc((lab+1)*t2->largeGrid->order, sizeof(double));
  
  for(pa=0; pa<Na; pa++) {
    zeta_a = zs_a[pa];
    s1 = 2.0*zeta_a*dAC;
    
    for(pb=0; pb<Nb; pb++) {
      zeta_b = zs_b[pb];
      s2 = 2.0*zeta_b*dBC;
      ps_Q.C = cs_a[pa] * cs_b[pb];
      /* derivatives: ps_Q.C gets an extra factor of zeta^n */
      for(n=0; n<shifta; n++)
	ps_Q.C *= zeta_a;
      for(n=0; n<shiftb; n++)
	ps_Q.C *= zeta_b;
      
      /* linear grid point mapping */
      grid = GCIntegrationTable_copy(t2->largeGrid);
      zeta_p = zeta_a + zeta_b;
      p = (zeta_a*dAC + zeta_b*dBC)/zeta_p;
      transformGrid_FM06(grid->order, grid->x, grid->w, zeta_p, p);
      
      /* tabulate r^(alpha+beta), bessel functions K(r), exp(-zeta*(d-r)**2), U(r) */
      for(n=0; n<grid->order; n++) {
	r = grid->x[n];
	/* formula for Q exponent modified as we are using weighted bessel functions 
	   cmp. Flores-Moreno et al., JCP 27 (2006), 1009--1019, eq. (24) */
	d1 = dAC - r;
	d2 = dBC - r;
	ps_Q.exponents[n] = - zeta_a*d1*d1 - zeta_b*d2*d2;
	/* simple screening: check for exp(t2->lnAccuracy) > 0.0 */
	if(r>dAC && r>dBC && ps_Q.exponents[n]<ps_Q.minExp) {
	  grid->end = n-1;
	  break;
	}
	else if(ps_Q.exponents[n]>=ps_Q.minExp) {
	  ps_Q.U[n] = evalECP(U, l, r);
	  weightedBesselFunction(t2->bessel, grid->order, laC, s1*r, K_a + n);
	  weightedBesselFunction(t2->bessel, grid->order, lbC, s2*r, K_b + n);
	  /* FIXME - when calculating r^(alpha+beta):
	     what is the convention for the exponent here - alpha+beta+2 or alpha+beta?
	     depends on whether the +2 has been included in the parameter set or not... */
	  rn[n] = 1.0;
	  for(l3=1; l3<=lab; l3++) {
	    rn[l3*grid->order + n] = r * rn[(l3-1)*grid->order + n];
	  }
	}
      }

      /* perform numerical integration with PS93 two-point sequence */      
      for(n=nFailed; n>=0; n--) {
	l1 = repeat.l1[n];
	l2 = repeat.l2[n];
	l3 = repeat.l3[n];
	ps_Q.K_a = K_a + l1*grid->order;
	ps_Q.K_b = K_b + l2*grid->order;
	ps_Q.rn = rn + l3*grid->order;
	
	if (0 != integrateGC_PSM92(&f_Q, grid)) {
	  /* error: numeric integration on large grid failed */
	  GCIntegrationTable_free(grid);
	  doubleArray_free(integrals);
	  goto endOfFunction;
	}
	else {
	  T[l1*Tinc2 + l2*Tinc1 + l3] += grid->I;
	  repeat.N--;
	}
      }

      GCIntegrationTable_free(grid);
    }
  } 
 endOfFunction:
  free(ps_Q.exponents);
  free(K_a);
  free(K_b);
  free(rn);
  free(ps_Q.U);
  return repeat;
}


/* combine the various contributions to the type2 integral */
doubleArray
calcGamma(Type2 *t2, doubleArray F, doubleArray UTab, ECP *U, const double dAC, const double dBC, 
	  const int la, const int shella, const int shellOffseta, doubleArray omegaA, doubleArray uspA, const int shifta,   
	  const int lb, const int shellb, const int shellOffsetb, doubleArray omegaB, doubleArray uspB, const int shiftb) {
  doubleArray T;  
  IntegrationFailure failed;
  int alpha, beta, lab;
  int lambda1, lambda2, l, m, c1, c2;
  double factor, tmp;
  int p, q;
  /* limits for sums over lambda */
  int ll1, ul1, ll2, ul2;
  int parity;
  int incT1, incT2;
  /* increments for omegaA,B arrays */
  const int incA1 = omegaA.dim[3];
  const int incA2 = omegaA.dim[2]*omegaA.dim[3];
  const int incB1 = omegaB.dim[3];
  const int incB2 = omegaB.dim[2]*omegaB.dim[3];
  /* increments for F array */
  const int Finc2 = F.dim[3]*F.dim[4];
  const int Finc3 = F.dim[2]*F.dim[3]*F.dim[4];
  /* integral storage */
  //doubleArray I = doubleArray_new(2, C_DIM(la), C_DIM(lb));
  doubleArray gamma = doubleArray_new(2, C_DIM(la), C_DIM(lb)); 
  const int inc = gamma.dim[2];

  /* calculate contributions due to angular and radial terms */
  for(l=0; l<U->L; l++) {
    T = doubleArray_new(3, la+l+1, lb+l+1, la+lb+1);
    /* 1a) calculate radial integrals */
    p = 0>shifta ? 0 : shifta;
    q = 0>shiftb ? 0 : shiftb;
    failed = calcT_fast_FM06(T, t2, l, UTab, la, lb, F,
			     &(F.array[p*Finc3+shella*Finc2]),
			     &(F.array[q*Finc3+shellb*Finc2]));
    if (0<failed.N)        
      failed = calcT_FM06(T, failed, t2, U, l, dAC, la, shella, shellOffseta, shifta, dBC, lb, shellb, shellOffsetb, shiftb);
    free(failed.l1);
    free(failed.l2);
    free(failed.l3);
    if (0<failed.N) {
      doubleArray_free(gamma);
      return gamma;
    } 
    incT1 = T.dim[3];	
    incT2 = T.dim[2]*T.dim[3];

    ul1 = la+l;
    ul2 = lb+l;
    /* 1b) link angular and radial terms */
    for(alpha=0; alpha<=la; alpha++) {
      /* conditions for the sums over lambda1,2 according to R. Flores-Moreno et al., J. Comput. Chem. 27 (2005), 1009.
         1.  |l-i-j-k| <= lambda <= l+i+j+k 
         2.  i+j+k+l-lambda must be even 
	 this implies:
         0 <= lambda
         i+j+k+l even(odd) => lambda even(odd) */          
      parity = (alpha+l)%2;
      ll1 = l-alpha;
      ll1 = (parity>ll1) ? parity : ll1;
      
      /* loop over cartesian components alpha_x,y,z */
      for(c1=0; c1<IJK_DIM(alpha); c1++) {
	p = C_INDEX(alpha,c1);
	// beta
	for(beta=0; beta<=lb; beta++) {
	  lab = alpha+beta;
	  parity = (beta+l)%2;
	  ll2 = l-beta;
	  ll2 = (parity>ll2) ? parity : ll2;
	  
	  // beta_{x,y,z}		    
	  for(c2=0; c2<IJK_DIM(beta); c2++) {
	    q = C_INDEX(beta,c2);
	    tmp = 0.0;
	    /* evaluate sums over lambda1,2 */
	    for(lambda1=ll1; lambda1<=ul1; lambda1+=2) {
	      for(lambda2=ll2; lambda2<=ul2; lambda2+=2) {
		factor = 0.0;
		for(m=0; m<2*l+1; m++) {
		  factor += omegaA.array[lambda1*incA2 + LM_INDEX(l,m)*incA1 + C_INDEX(alpha,c1)] 
		          * omegaB.array[lambda2*incB2 + LM_INDEX(l,m)*incB1 + C_INDEX(beta,c2)];		  
		}
		tmp += factor*T.array[lambda1*incT2 + lambda2*incT1 + lab];	   
	      }
	    }
	    gamma.array[p*inc+q] += tmp;
	  }
	}
      }
    }
    doubleArray_free(T);  
  }
  
  /* 2. transform polynomials into local coordinate system */
  return calcPolynomials(gamma, t2->ijk, t2->ijkIndex, t2->ijkDim,
			 t2->fac, t2->accuracy, 16.0*M_PI*M_PI, la, uspA, lb, uspB);
}

