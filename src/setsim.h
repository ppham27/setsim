#ifndef SETSIM_H_
#define SETSIM_H_

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>

typedef struct {
  SEXP h;
  SEXP env;
  SEXP epsilon;
  SEXP input;
  SEXP vz;
  SEXP cstar;
  SEXP fcall;
} hFunction;

double * calculateVZMatrix(int mv, int nv, double *v,
                           int mz, int nz, double *z);

int nroots(int nh, double *h,
           int ne, double *e,
           int *position,
           double *range);

double bisect(hFunction *f,
              double a,
              double b);

SEXP boundaryFixed(SEXP h, SEXP hEnv,
                   SEXP v,
                   SEXP z,
                   SEXP zt,
                   SEXP critValue,
                   SEXP targetValue,
                   SEXP input,
                   SEXP mle,
                   SEXP solType);


#endif /* SETSIM_H_ */
