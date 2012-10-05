
#ifndef FACTORIAL_H
#define FACTORIAL_H 1

/* calculate n! ineffectively */
double * factorials(const int n);

/* calculate n!! ineffectively (even and odd numbers n) */
double * double_factorials(const int n);

/* calculate the binomial coefficient n over k */
double n_over_k(const int n, const int k, const double *fac);

#endif
