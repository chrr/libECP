/* Copyright (c) 2012, Christoph Reimann */

#ifndef GC_INTEGRATORS_H
#define GC_INTEGRATORS_H

typedef double (*IntegrandFunction)(double x, int index, const void *params);

typedef struct {
  /* callback function f takes 3 arguments:
     - abscissa x
     - index of x in the GC grid
     - arbitrary parameters */
  IntegrandFunction f;
  void *params;
} Integrand;

typedef struct {
  int n; 
  /* number of points */
  int order; 
  /* abscissae */
  double *x; 
  /* weights */
  double *w; 
  double tol;
  double I;
  /* grid dimensions in case of screening */
  int start, end;
} GCIntegrationTable;

GCIntegrationTable * GCIntegrationTable_copy(const GCIntegrationTable *O);
void GCIntegrationTable_free(GCIntegrationTable *t);

/* Gauss-Chebyshev quadrature of the second kind according to Perez-Jordan, San-Fabian, Moscardo (PSM)
   J.M. Perez-Jorda et al., Comput. Phys. Comm. 70 (1992), 271--284. */
GCIntegrationTable * initGridGC_PSM92(int maxPoints, const double tol);
int integrateGC_PSM92(const Integrand *F, GCIntegrationTable *t);

/* Gauss-Chebyshev quadrature of the second kind according to Perez-Jordan and San-Fabian (PS)
   J.M. Perez-Jorda et al., Comput. Phys. Comm. 77 (1993), 46--56.
   
   features compared to PSM92:
   + total number of points  required for computing an integral increases more moderately
   + more conservative error estimate
   - generally somewhat less efficient (regarding the total number of points required),
     but depends strongly on the integrand */
GCIntegrationTable * initGridGC_PS93(int maxPoints, const double tol);
int integrateGC_PS93(const Integrand *F, GCIntegrationTable *t);

/* transform abscissas from (-1,+1) to (0,\infty)
   M. Krack, A.M. K"oster, J. Chem. Phys. 108 (1998), 3226--323.4  */
void transformGrid_KK(int n, double *x, double *w);

/* transfrom abscissas from (-1,+1) to (rmin,rmax) using a linear mapping
   R. Flores-Moreno, R.J. Alvarez-Mendez, A. Vela, A.M. K"oster, J. Comput. Chem. 27 (2006), 1009--1019. */
void transformGrid_FM06(int n, double *x, double *w, double zeta_P, double P);

#endif
