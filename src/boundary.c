#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "boundary.h"

SEXP boundaryFixed(SEXP h, SEXP hEnv,
                   SEXP v,
                   SEXP z,
                   SEXP zt,
                   SEXP critValue,
                   SEXP targetValue,
                   SEXP input,
                   SEXP mle,
                   SEXP solType) {
  int mv = INTEGER(getAttrib(v, R_DimSymbol))[0];
  int nv = INTEGER(getAttrib(v, R_DimSymbol))[1];
  int mz = INTEGER(getAttrib(z, R_DimSymbol))[0];
  int nz = INTEGER(getAttrib(z, R_DimSymbol))[1];
  int mzt = INTEGER(getAttrib(zt, R_DimSymbol))[0];
  int nzt = INTEGER(getAttrib(zt, R_DimSymbol))[1];
  int mSolType = INTEGER(getAttrib(solType, R_DimSymbol))[0];
  int nSolType = INTEGER(getAttrib(solType, R_DimSymbol))[1];

  double *vzMatrix = calculateVZMatrix(mv, nv, REAL(v),
                                       mz, nz, REAL(z));
  int mvz = mv;
  int nvz = mz;

  /* initialize function */
  hFunction *hFunc;
  hFunc = (hFunction *) R_alloc(1, sizeof(hFunction));
  PROTECT(hFunc -> h = h);
  PROTECT(hFunc -> env = hEnv);
  PROTECT(hFunc -> epsilon = allocVector(REALSXP, 1));
  REAL(hFunc -> epsilon)[0] = 0;
  PROTECT(hFunc -> input = input);
  PROTECT(hFunc -> vz = allocVector(REALSXP, mvz));
  PROTECT(hFunc -> cstar = allocVector(REALSXP, 1));
  REAL(hFunc -> cstar)[0] = 0;
  PROTECT(hFunc -> fcall = lang5(hFunc -> h,
                                 hFunc -> epsilon,
                                 hFunc -> input,
                                 hFunc -> vz,
                                 hFunc -> cstar));

  /* initiate values for nroots function */
  double *epsilon;
  double *hValues;
  int *roots_position;
  double *roots_range;  
  int nEpsilon = 51;
  int nhValues = nEpsilon;
  int nRoots_range = 4;
  epsilon = (double *) R_alloc(nEpsilon, sizeof(double));
  hValues = (double *) R_alloc(nhValues, sizeof(double));
  roots_position = (int *) R_alloc(nhValues, sizeof(int));
  roots_range = (double *) R_alloc(nhValues*nRoots_range, sizeof(double));
  
  double a, zSquared;
  int i, j, k;
  double *out;
  double *wald;
  int nOut = 4+mvz;
  int mWald = mz*2;
  int nWald = 4+mvz;

  out = (double *) R_alloc(100000, sizeof(double));
  wald = (double *) R_alloc(mWald*nWald, sizeof(double));
  int mOut = 0;
  for (i = 0; i < mz; i++) {
    int roots_n;

    memcpy(REAL(hFunc -> vz), vzMatrix + i*mvz, mvz*sizeof(double));
    for (j = 0; j < LENGTH(critValue); j++) {
      REAL(hFunc -> cstar)[0] = REAL(critValue)[j];
      /* calculate a */
      zSquared = 0;
      for (k = 0; k < mzt; k++) {
        zSquared += REAL(zt)[i*mzt + k]*REAL(zt)[i*mzt + k];
      }
      a = sqrt(REAL(hFunc -> cstar)[0]/zSquared);
      
      /* calculate epsilon and hValues in preparation of n roots */
      epsilon[0] = -4*a;
      REAL(hFunc -> epsilon)[0] = epsilon[0];
      hValues[0] = REAL(eval(hFunc -> fcall,
                             hFunc -> env))[0];
      for (k = 1; k < nEpsilon; k++) {
        epsilon[k] = epsilon[k-1] + 8*a/(nEpsilon-1);
        REAL(hFunc -> epsilon)[0] = epsilon[k];
        hValues[k] = REAL(eval(hFunc -> fcall,
                               hFunc -> env))[0];
      }
      roots_n = nroots(nhValues, hValues,
                       nEpsilon, epsilon,
                       roots_position,
                       roots_range);

      if (roots_n == 0) {
        REAL(solType)[j + mSolType] += 1;
      } else {
        double sol;

        for (int r = 0; r < roots_n; r++) {
          sol = bisect(hFunc,
                       roots_range[r*nRoots_range],
                       roots_range[r*nRoots_range + 1]);
          out[mOut*nOut] = i + 1;
          out[mOut*nOut + 1] = REAL(targetValue)[j];
          out[mOut*nOut + 2] = roots_n;
          out[mOut*nOut + 3] = sol;
          for (int k = 0; k < mvz; k++) {
            out[mOut*nOut + 4 + k] = REAL(mle)[k] + sol*REAL(hFunc -> vz)[k];
          }
          mOut++;
        }
        switch (roots_n) {
        case 1:
          REAL(solType)[j + 2*mSolType] += 1;
          break;
        case 2:
          REAL(solType)[j + 3*mSolType] += 1;
          break;
        default:
          REAL(solType)[j + 4*mSolType] += 1;
          break;
        }
      }
    }
    double sol_w = sqrt(3.841459/zSquared);
    wald[2*i*nWald] = i + 1;
    wald[2*i*nWald + 1] = REAL(targetValue)[j-1];
    wald[2*i*nWald + 2] = roots_n;
    wald[2*i*nWald + 3] = -sol_w;
    for (int k = 0; k < mvz; k++) {
      wald[2*i*nWald + 4 + k] = REAL(mle)[k] - sol_w*REAL(hFunc -> vz)[k];
    }
    wald[(2*i+1)*nWald] = i + 1;
    wald[(2*i+1)*nWald + 1] = REAL(targetValue)[j-1];
    wald[(2*i+1)*nWald + 2] = roots_n;
    wald[(2*i+1)*nWald + 3] = sol_w;
    for (int k = 0; k < mvz; k++) {
      wald[(2*i+1)*nWald + 4 + k] = REAL(mle)[k] + sol_w*REAL(hFunc -> vz)[k];
    }
    if (i % 200 == 0) {
      Rprintf("Generated %d rays.\n", i+1);
    }
  }
  Rprintf("Done!\n");
  SEXP output;
  SEXP out_sexp;
  SEXP wald_sexp;  
  PROTECT(output = allocVector(VECSXP, 2));
  PROTECT(out_sexp = allocVector(REALSXP, mOut*nOut));
  PROTECT(wald_sexp = allocVector(REALSXP, mWald*nWald));
  memcpy(REAL(out_sexp), out, sizeof(double)*mOut*nOut);
  memcpy(REAL(wald_sexp), wald, sizeof(double)*mWald*nWald);
  SET_VECTOR_ELT(output, 0, out_sexp);
  SET_VECTOR_ELT(output, 1, wald_sexp);
  UNPROTECT(10);
  return output;
}

double * calculateVZMatrix(int mv, int nv, double *v,
                           int mz, int nz, double *z) {
  double *c;
  c = (double *) R_alloc(mv*mz,sizeof(double));

  double alpha = 1;
  double beta = 0;
  F77_CALL(dgemm)("n","t", &mv, &mz, &nv, &alpha,
                  v, &mv,
                  z, &mz,
                  &beta,
                  c, &mv);
  return c;
}

int nroots(int nh, double *h,
           int ne, double *e,
           int *position,
           double *range) {
  int i;
  int n = 0;
  for (i = 1; i < nh; i++) {
    if (h[i] == 0 || (h[i-1] < 0 && h[i] > 0) || (h[i-1] > 0 && h[i] < 0)) {
      position[n] = i;
      range[n*4] = e[i-1];
      range[n*4 + 1] = e[i];
      range[n*4 + 2] = h[i-1];
      range[n*4 + 3] = h[i];
      n++;
    } 
  }
  return n;
}

double bisect(hFunction *f,
              double a,
              double b) {
  double res = 0;
  double diff = 10;
  double tol = 0.00001;
  double oldP = (a+b)/2.0;
  REAL(f -> epsilon)[0] = a;
  double f_a = REAL(eval(f -> fcall,
                         f -> env))[0];
  if (f_a==0) { return a; }
  REAL(f -> epsilon)[0] = b;
  double f_b = REAL(eval(f -> fcall,
                         f -> env))[0];
  if (f_b==0) { return b; }
  double f_old_p;
  double f_new_p;
  double newP;
  int i=1;
  while (diff >= tol) {
    REAL(f -> epsilon)[0] = oldP;
    f_old_p = REAL(eval(f -> fcall,
                        f -> env))[0];
    REAL(f -> epsilon)[0] = a;    
    f_a = REAL(eval(f -> fcall,
                    f -> env))[0];
    if (f_a * f_old_p > 0) {
      a = oldP;
    } else {
      b = oldP;
    }
    newP = (a+b)/2.0;
    REAL(f -> epsilon)[0] = newP;
    f_new_p = REAL(eval(f -> fcall,
                        f -> env))[0];
    diff = fabs(f_new_p - f_old_p);
    oldP = newP;
    res = oldP;
    i++;
  }
  return res;
}
