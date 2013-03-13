/* Copyright (c) 2012, Christoph Reimann */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "gc_integrators.h"
#include "dimensions.h"
#include "transformations.h"


int
inefficient_factorial(int new_n, double **fac, int increment) {
  int i;
  double *f;

  assert(new_n >= 0);
  assert(1==increment || 2==increment);

  *fac = calloc(new_n+1, sizeof(double));
  f = *fac;
  if (0<new_n)
    f[0] = 1.0;
  if (1<new_n)
    f[1] = 1.0;
  for(i=2; i<=new_n; i++) {
    f[i] = f[i-increment]*i;
    if (isinf(f[i]))
      abort();
  }
  
  return new_n;
}


/* calculate n! for 0..n */
double *
factorials(const int n) {
  double *f = NULL;

  inefficient_factorial(n, &f, 1);

  return f;
}


/* calculate n!! for 0..n */
double *
double_factorials(const int n) {
  double *f = NULL;

  inefficient_factorial(n, &f, 2);
      
  return f;
}


/* function for calculating the binomial coefficient n over k for 
   0 <= k <= n (zero otherwise) */
double
n_over_k (const int n, const int k, const double *fac) {
  double result = 0.0;
    
  if (k>=0 && k<=n) {
    result = fac[n]/(fac[n-k]*fac[k]);
  }
  
  return result;
}


/* convert cartesian to spherical coordinates */
void
sphericalCoordinates (double *cartesian, double *spherical) {
  double r, theta, phi;
  const double accuracy = 1.0E-14;
  const double x = cartesian[0];
  const double y = cartesian[1];
  const double z = cartesian[2];

  r = sqrt(x*x+y*y+z*z);
  if (r < accuracy)
    theta = 0.0;
  else
    theta = acos(z/r);

  if (fabs(x) < accuracy) {
    if (fabs(y) < accuracy) 
      phi = 0.0;
    else if (y < 0.0)
      phi = 1.5*M_PI;
    else 
      phi = 0.5*M_PI;
  }
  else {
    if (x > 0.0) 
      phi = atan(y/x);
    else 
      phi = atan(y/x) + M_PI;
  }

  spherical[0] = r;
  spherical[1] = theta;
  spherical[2] = phi;
}


double
distance(double *A, double *B) {
  double x,y,z;
  x = A[0]-B[0];
  y = A[1]-B[1];
  z = A[2]-B[2];
  return sqrt(x*x+y*y+z*z);
}


void
distanceVector(double *BA, double *A, double *B) {
  BA[0] = B[0] - A[0];
  BA[1] = B[1] - A[1];
  BA[2] = B[2] - A[2];
}


/* calculate atomic radius for primitive gaussian function
   parameters:
   c: contraction coefficient
   zeta: exponent
   l: radial power
   threshold: cutoff value for radius */
double
gtoRadius (const double c, const double zeta, const int l, const double cutoff) {
  int i;
  double delta, dg, guess, zr, r = 0.0;
  const int N = 40;
  int converged = 0;
  const double t = log(fabs(c)/fabs(cutoff));
  
  /* analytical expression for l = 0 */
  dg = t/zeta;
  r = (dg>cutoff) ? dg : cutoff;

  /* l>0: find outer r by Newton-Raphson procedure */
  if (0!=l) {
    /* for large radial powers: correct r if necessary */
    if (0<l) {
      dg = sqrt(0.5*l/fabs(c));
      guess = (dg>cutoff) ? dg : cutoff;
      if (guess > r) {
	r = 0.5*(r+guess);
      }
    }
    for(i=0; i<N; i++) {
      zr = zeta*r;
      guess = t + l*log(r) - zr*r;
      delta = guess/(l/r - 2*zr);
      dg = r-delta;
      r = (dg>cutoff) ? dg : cutoff;

      if (fabs(delta) < cutoff) {
	r *= r;
	converged = 1;
	break;
      }
    }
    /* possible error: no outer atomic radius found for primitive gaussian */
    assert(0 != converged);
  }
  
  return r;
}


/* calculate radial cutoff for a gaussian shell */
double
shellRadius(int depth, int am, double *d, double *a, const double zero) {
  int i;
  double zeta = a[0];
  double c = fabs(d[0]);
  double radius = 0.0;

  for(i=1; i<depth; i++) {
    if (a[i]<zeta && d[i]!=0.0) {
      zeta = a[i];
      c = fabs(d[i]);
    }	    
  }
  radius = gtoRadius(c, zeta, am, zero);
  return sqrt(radius);
}


/* screening by potential cutoff
   returns: index (elements U[0..index] are to be retained)
   -1: potential is screened completely */
int
potentialScreening(GCIntegrationTable *grid, double *U, double zero) {
  int i, cutoff_index = -1;
  /* U[i <= cutoff_index] is to be treated non-zero */
  for(i=grid->end; i>=grid->start; i--) {
    if(fabs(U[i]) > zero) {
      cutoff_index = i;
      break;
    }
  }
    
  return cutoff_index;
}


/* evaluate unitary sphere polynomial terms x^i y^j z^k */
doubleArray
unitarySpherePolynomials(const int l, double *r) {
  int i, j, k;
  double x, y, z;
  const int dim = l+1;
  const int inc1 = dim;
  const int inc2 = dim*dim;
  doubleArray table = doubleArray_new(3, dim, dim, dim);

  for(i=0; i<=l; i++) {
    if (0 == i)
      x = 1.0;
    else
      x = -x*r[0];
    for(j=0; j<=l-i; j++) {
      if (0 == j) 
	y = 1.0;
      else
	y = -y*r[1];
      for(k=0; k<=l-i-j; k++) {
	if (0 == k) 
	  z = 1.0;
	else
	  z = -z*r[2];
	table.array[i*inc2 + j*inc1 + k] = x*y*z;
      }
    }
  }
  return table;
}


doubleArray
calcPolynomials(doubleArray gamma, const int *ijk, const int *ijkIndex, const int ijkDim, 
		double *fac, const double tol, const double N, 
		const int la, doubleArray uspA, 
		const int lb, doubleArray uspB) {
  int c1, c2, beta;
  int ax, ay, az, bx, by, bz;
  int alpha_x, alpha_y, alpha_z, beta_x, beta_y, beta_z;
  int dax, day, daz, dbx, dby, dbz;
  double bin_ax, bin_ay, bin_az, bin_bx, bin_by, bin_bz;
  int p, q, r;
  double factor;
  /* array increments */
  const int incG = gamma.dim[2];
  const int incA1 = uspA.dim[3];
  const int incA2 = uspA.dim[2]*uspA.dim[3];
  const int incB1 = uspB.dim[3];
  const int incB2 = uspB.dim[2]*uspB.dim[3];
  doubleArray J = doubleArray_new(2, IJK_DIM(la), C_DIM(lb));
  doubleArray I = doubleArray_new(2, IJK_DIM(la), IJK_DIM(lb));
  const int incJ = J.dim[2];
  const int incI = I.dim[2];

  /* a) sums over alpha_x,y,z */
  for(c1=0; c1<IJK_DIM(la); c1++) {
    p = CIJK_INDEX(la,c1);
    ax = ijk[p];
    ay = ijk[p+1];
    az = ijk[p+2];
    for(alpha_x=0; alpha_x<=ax; alpha_x++) {
      bin_ax = n_over_k(ax, alpha_x, fac);
      dax = ax - alpha_x;
      for(alpha_y=0; alpha_y<=ay; alpha_y++) {
	bin_ay = bin_ax * n_over_k(ay, alpha_y, fac);
	day = ay - alpha_y;
	for(alpha_z=0; alpha_z<=az; alpha_z++) {
	  bin_az = bin_ay * n_over_k(az, alpha_z, fac);
	  daz = az - alpha_z;

	  factor = bin_az * uspA.array[dax*incA2 + day*incA1 + daz];
	  if (fabs(factor) <= tol)
	    continue;

	  r = ijkIndex[alpha_x*ijkDim*ijkDim + alpha_y*ijkDim + alpha_z];
	  for(beta=0; beta<=lb; beta++) {
	    for(c2=0; c2<IJK_DIM(beta); c2++) {
	      q = C_INDEX(beta,c2);
	      J.array[c1*incJ+q] += factor*gamma.array[r*incG+q]; 
	    }
	  }
	}
      }
    }
  }

  /* b) sums over beta_x,y,z */
  for(c2=0; c2<IJK_DIM(lb); c2++) {
    q = CIJK_INDEX(lb,c2);
    bx = ijk[q];
    by = ijk[q+1];
    bz = ijk[q+2];
    for(beta_x=0; beta_x<=bx; beta_x++) {
      bin_bx = n_over_k(bx, beta_x, fac);
      dbx = bx - beta_x;
      for(beta_y=0; beta_y<=by; beta_y++) {
	bin_by = bin_bx * n_over_k(by, beta_y, fac);
	dby = by - beta_y;
	for(beta_z=0; beta_z<=bz; beta_z++) {
	  bin_bz = bin_by * n_over_k(bz, beta_z, fac);
	  dbz = bz - beta_z;

	  factor = bin_bz*uspB.array[dbx*incB2 + dby*incB1 + dbz];
	  if (fabs(factor) <= tol)
	    continue;

	  factor *= N;
	  r = ijkIndex[beta_x*ijkDim*ijkDim + beta_y*ijkDim + beta_z];
	  for(c1=0; c1<IJK_DIM(la); c1++) {
	    I.array[c1*incI+c2] += factor*J.array[c1*incJ+r];
	  }
	}
      }
    }
  }

  doubleArray_free(gamma);
  doubleArray_free(J);
  return I;
}

