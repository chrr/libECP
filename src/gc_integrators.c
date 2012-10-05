
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "gc_integrators.h"

GCIntegrationTable *
GCIntegrationTable_copy(const GCIntegrationTable *O) {
  int i;
  GCIntegrationTable *C = calloc(1, sizeof(GCIntegrationTable));
  C->order = O->order;
  C->n = O->n;
  C->x = calloc(C->order, sizeof(double));
  C->w = calloc(C->order, sizeof(double));
  C->tol = O->tol;
  C->I = O->I;
  C->start = O->start;
  C->end = O->end;
  for(i=0; i<C->order; i++) {
    C->x[i] = O->x[i];
    C->w[i] = O->w[i];
  }
  return C;
}

void 
GCIntegrationTable_free(GCIntegrationTable *t) {
  free(t->x);
  free(t->w);
  free(t);
}


/* Gauss-Chebyshev quadrature of the second kind according to Perez-Jordan, San-Fabian, Moscardo (PSM)
   J.M. Perez-Jorda et al., Comput. Phys. Comm. 70 (1992), 271--284. */
int
integrateGC_PSM92(const Integrand *F, GCIntegrationTable *t) {
  /* check the maximal number: log(n-1)/log(2) must be an integer */
  int order = pow(2, floor(log(t->n+1)/log(2)))-1;
  double *x = t->x;
  double *w = t->w;
  int runs = (int)floor(log(order)/log(2));
  /* offset in x and w arrays */
  int offset = (int)pow(2, runs);
  const int M = (order-1)/2;
  /* initialize integral */
  double I = w[M]*F->f(x[M], M, F->params);
  int n = 1;
  double N = n+1.0;
  double e, T, q, p = I;
  int idx, i, cnt;

  while(n<=M) {
    q = 2*p;
    p = 2*I;
    offset /= 2;
    cnt = 0;
    for(i=1; i<=n; i+=2) {
      idx = i*offset-1;
      T = 0.0;
      if(idx >= t->start) {
	T += w[idx]*F->f(x[idx], idx, F->params);
	cnt++;
      }
      if(order-idx-1 <= t->end) {
	T += w[order-idx-1]*F->f(x[order-idx-1], order-idx-1, F->params);
	cnt++;
      }
      I += T;
    }
    n = 2*n+1;
    N = n+1.0;
    e = I-p;
    if(0==cnt) 
      continue;

    /* check for convergence */
    if(16*e*e <= 3*N*fabs(I-q)*t->tol) {
      t->I = 16*I/(3*N);
      return 0;
    }      
  }
  return 1;
}


GCIntegrationTable *
initGridGC_PSM92(const int maxPoints, const double tol) {
  /* determine the max. order: log(n-1)/log(2) must be an integer */
  int order = pow(2, floor(log(maxPoints+1)/log(2)))-1;
  int runs = (int)floor(log(order)/log(2));
  /* offset in x and w arrays */
  int offset = (int)pow(2, runs);
  int M = (order-1)/2;
  int n = 1;
  double N = n+1.0;
  int i;
  int idx;
  double S0 = 1.0;
  double C0 = 0.0;
  double S1, C1, s, c, t;
  GCIntegrationTable *T = calloc(1, sizeof(GCIntegrationTable));
  double *x = T->x = calloc(order, sizeof(double));
  double *w = T->w = calloc(order, sizeof(double));
  T->n = T->order = order;
  T->start = 0;
  T->end = order-1;
  T->tol = tol;

  /* x(1,1)=0, w(1,1)~sin(pi/2)=1 reside in the middle of the interval */
  x[M] = 0.0;
  w[M] = 1.0;

  while(n <= M) {
    C1 = C0;
    S1 = S0;
    C0 = sqrt((1+C1)/2);
    S0 = S1/(2*C0);
    s = S0;
    c = C0;

    /* make use of the relations
          x(2i, 2n+1) == x(i,n)
          w(2i, 2n+1) == 0.5 w(i,n) */
    offset /= 2;
    for(i=1; i<=n; i+=2) {
      t = 1 + 2/(3*M_PI) * (3+2*s*s)*s*c - i/N;
      idx = i*offset-1;
      /* +x */
      x[order-idx-1] = t;
      /* -x */
      x[idx] = -t;
      w[order-idx-1] = w[idx] = s*s*s*s;
      t = s;
      s = s*C1 + c*S1;
      c = c*C1 - t*S1;
    }
    n = 2*n+1;
    N = n+1.0;
  }

  return T;
}


/* Gauss-Chebyshev quadrature of the second kind according to Perez-Jordan and San-Fabian (PS)
   J.M. Perez-Jorda et al., Comput. Phys. Comm. 77 (1993), 46--56.
   
   features compared to PSM92:
   + total number of points  required for computing an integral increases more moderately
   + more conservative error estimate
   - generally somewhat less efficient (regarding the total number of points required),
     but depends strongly on the integrand */
int
integrateGC_PS93(const Integrand *F, GCIntegrationTable *t) {
  int runs = (int)floor(log(t->n)/log(2));
  /* offset in x and w arrays */
  int offset = (int)pow(2, runs);
  int order = 3*offset-1;
  double *x = t->x;
  double *w = t->w;
  int n = 3;
  /* j is a switch: 1 => two-point formula
                    0 => one-point formula */
  int j = 0;
  int i, cnt, idx = offset-1;
  double p, q, I, err, T;
  const void *ps = F->params;

  assert(order>=3);

  /* initialize integral: I[1,1], I[1,2] and I[2,2] */
  p = w[(order-1)/2] * F->f(0.0, (order-1)/2, ps);
  q = w[idx] * F->f(x[idx], idx, ps) + w[order-offset] * F->f(x[order-offset], order-offset, ps);
  I = p + q;

  offset /= 2;
  
  while((2*n*(1-j) + j*4*n/3 - 1) <= order) {
    j = 1-j;
    if (0==j)
      offset /= 2;
    cnt = 0;
    for(i=1; i<n; i+=2) {
      if (3*((i+2*j)/3) >= i+j) {
	idx = i*offset-1;
	T = 0.0;
	if(idx >= t->start) {
	  T += w[idx] * F->f(x[idx], idx, ps);
	  cnt++;
	}
	if (order-idx-1 <= t->end) {
	  T += w[order-idx-1] * F->f(x[order-idx-1], order-idx-1, ps);
	  cnt++;
	}
	I += T;
      }
    }
    i = n;
    n *= (1+j);
    p += (1-j)*(I-q);
    if(0<cnt) 
      err = 16*fabs((1-j)*(q-3*p/2) + j*(I-2*q)) / (3*n);
    q = (1-j)*q + j*I;
    if(0==cnt) 
      continue;    

    /* check for convergence */
    if (err < t->tol) {
      t->I = 16*q/(3*n);
      return 0;
    }      
  }
  return 1;
}


GCIntegrationTable *
initGridGC_PS93(const int maxPoints, const double tol) {
  int runs = (int)floor(log(maxPoints)/log(2));
  int offset = (int)pow(2, runs);
  int i, idx, n, order;
  double C0, S0, S1, C1, s, s2, c, t;
  double *x, *w;
  GCIntegrationTable *T = calloc(1, sizeof(GCIntegrationTable));
  
  /* include 2m+1 abscissas + x=0.0 */
  T->order = order = 3*offset - 1;
  assert(order>=3);
  T->start = 0;
  T->end = order-1;
  /* values for runs, offset, order will be recalculated in integrateGC_PS93 */
  T->n = maxPoints;
  T->tol = tol;
  x = T->x = calloc(order, sizeof(double));
  w = T->w = calloc(order, sizeof(double));

  n = 3;
  C0 = sin(M_PI/3);
  S0 = 0.5;
  C1 = S0;
  S1 = C0;
  /* I[1,2] and I[2,2] */
  c = cos(M_PI/3);
  s = C0;
  s2 = s*s;
  
  /* x(1,1)=0, w(1,1)~sin(pi/2)=1 reside in the middle of the interval */
  x[order/2] = 0.0;
  w[order/2] = 1.0;
  t = (n-2.0)/n + 2/M_PI*(1+2*s2/3)*c*s;
  x[offset-1] = -t;
  x[order-offset] = t;
  w[order-offset] = w[offset-1] = s2*s2;
        
  while ((4*n/3 - 1) <= order) {
    c = C0;
    s = S0;
    offset /= 2;
    
    for(i=1; i<n; i+=2) {
      s2 = s*s;
      idx = i*offset-1;
      t = 1 + 2/(3*M_PI)*s*c*(3+2*s2) - ((double)i)/n;
      x[idx] = -t;
      x[order-idx-1] = t;
      w[order-idx-1] = w[idx] = s2*s2;
      
      t = s;
      s = s*C1 + c*S1;
      c = c*C1 - t*S1;
    }
    n *= 2;
    C1 = C0;
    S1 = S0;
    C0 = sqrt((1+C0)/2);
    S0 = S0/(2*C0);
  }

  return T;
}


/* As the abscissas generated according to the GC scheme
   lie in the interval [-1,+1], they have to be transformed to the
   interval [0,+\infty).
   For this purpose Becke originally proposed the transformation
   r_A = R_A (1+x_A)/(1-x_A)  <=>  x_A = (r_A-R_A)/(r_A+R_A)
   where R_A is the 'bonding region'.
   Treutler and Ahlrichs (TA) suggested mappings of logarithmic type
   that projected more points into the chemically relevant bonding
   region and less points at large distances.
   Krack and K"oster (KK) modified TA's logarithmic transformation by
   omitting the atomic scaling parameter R_A, thus obtaining the
   simple and parameter-free transformation
   r_A = 1/ln2 ln(2/(1-x_A))  <=>  x_A = 1-2^(1-r_A) 
   
   M. Krack, A.M. K"oster, J. Chem. Phys. 108 (1998), 3226--3234 */
void
transformGrid_KK(int n, double *x, double *w) {
  int i;
  double ln2 = log(2.0);
  double x_i, w_i;

  for (i=0; i<n; i++) {
    x_i = 1.0 - log(1.0-x[i])/ln2;
    w_i = w[i]/(ln2*(1.0-x[i]));
    x[i] = x_i;
    w[i] = w_i;
  }
}


void 
transformGrid_FM06(int n, double *x, double *w, double zeta_P, double P) {
  double sigma = 1.0/sqrt(zeta_P);
  double t = P-7.0*sigma;
  double rmin = (t>0.0) ? t : 0.0;
  double rmax = P+9.0*sigma; 
  double i1 = 0.5*(rmax-rmin);
  double i2 = 0.5*(rmax+rmin);
  int i;

  for(i=0; i<n; i++) {
    x[i] = i1*x[i] + i2;
    /* weight factor is simply scaled */
    w[i] *= i1;
  }
}
