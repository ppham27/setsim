#include "setsim.h"

SEXP independentFixed(SEXP h, SEXP hEnv,
                      SEXP hLik, SEXP hLikEnv,
                      SEXP v,
                      SEXP z,
                      SEXP zt,
                      SEXP cstar,
                      SEXP targetValue,
                      SEXP input,
                      SEXP est,
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

  /* initialize functions */
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

  hLikFunction *hLikFunc;
  hLikFunc = (hLikFunction *) R_alloc(1, sizeof(hLikFunction));
  PROTECT(hLikFunc -> hLik = hLik);
  PROTECT(hLikFunc -> env = hLikEnv);
  PROTECT(hLikFunc -> y = VECTOR_ELT(input, getListIndex(input,"y")));
  PROTECT(hLikFunc -> X = VECTOR_ELT(input, getListIndex(input,"X")));
  PROTECT(hLikFunc -> est = allocVector(REALSXP, LENGTH(est)));
  memcpy(REAL(hLikFunc -> est), REAL(est), LENGTH(est)*sizeof(double));
  PROTECT(hLikFunc -> fcall = lang4(hLikFunc -> hLik,
                                    hLikFunc -> y,
                                    hLikFunc -> X,
                                    hLikFunc -> est));
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

  /* allocate output */
  double *out;
  double *wald;
  size_t outSize = 1000000;
  int nOut = 5+mvz;
  int mWald = mz*2;
  int nWald = 4+mvz;
  out = (double *) Calloc(outSize, double);
  wald = (double *) Calloc(mWald*nWald, double);
  int mOut = 0;

  double sol_w = 1;
  int i, j;
  int roots_n;
  double sol;
  for (i = 0; i < mz; i++) {
    /* start calculations */
    memcpy(REAL(hFunc -> vz), vzMatrix + i*mvz, mvz*sizeof(double));
    REAL(hFunc -> cstar)[0] = REAL(cstar)[i];
    /* calculate epsilon and hValues */
    epsilon[0] = -4.0;
    REAL(hFunc -> epsilon)[0] = epsilon[0];
    hValues[0] = REAL(eval(hFunc -> fcall,
                           hFunc -> env))[0];
    for (j = 1; j < nEpsilon; j++) {
      epsilon[j] = epsilon[j-1] + 8.0/(nEpsilon-1);
      REAL(hFunc -> epsilon)[0] = epsilon[j];
      hValues[j] = REAL(eval(hFunc -> fcall,
                             hFunc -> env))[0];
    }
    roots_n = nroots(nhValues, hValues,
                     nEpsilon, epsilon,
                     roots_position,
                     roots_range);
    if (roots_n == 0) {
      REAL(solType)[0] += 1;
      out[mOut*nOut] = i + 1;
      out[mOut*nOut + 1] = REAL(targetValue)[i];
      out[mOut*nOut + 2] = roots_n;
      for (j = 3; j < nOut; j++) { out[mOut*nOut + j] = R_NaReal; }
      mOut++;
    } else {
      for (int r = 0; r < roots_n; r++) {
        sol = bisect(hFunc,
                     roots_range[r*nRoots_range],
                     roots_range[r*nRoots_range + 1]);
        if (mOut*nOut + nOut >= outSize) {
          outSize *= 2;
          out = (double *) Realloc(out, outSize, double);
        }
        out[mOut*nOut] = i + 1;
        out[mOut*nOut + 1] = REAL(targetValue)[i];
        out[mOut*nOut + 2] = roots_n;
        out[mOut*nOut + 3] = sol;
        for (j = 0; j < mvz; j++) {
          REAL(hLikFunc -> est)[j] = out[mOut*nOut + 4 + j] = REAL(est)[j] + sol*REAL(hFunc -> vz)[j]; 
        }
        out[mOut*nOut + nOut - 1] = REAL(eval(hLikFunc -> fcall, hLikFunc -> env))[0];
        mOut++;
      }
      switch (roots_n) {
      case 1:
        REAL(solType)[1] += 1;
        break;
      case 2:
        REAL(solType)[2] += 1;
        break;
      default:
        REAL(solType)[3] += 1;
        break;
      }
    }
    wald[2*i*nWald] = i + 1;
    wald[2*i*nWald + 1] = REAL(cstar)[i];
    wald[2*i*nWald + 2] = roots_n;
    wald[2*i*nWald + 3] = -sol_w;
    for (j = 0; j < mvz; j++) {
      wald[2*i*nWald + 4 + j] = REAL(est)[j] - sol_w*REAL(hFunc -> vz)[j];
    }
    wald[(2*i+1)*nWald] = i + 1;
    wald[(2*i+1)*nWald + 1] = REAL(cstar)[i];
    wald[(2*i+1)*nWald + 2] = roots_n;
    wald[(2*i+1)*nWald + 3] = sol_w;
    for (j = 0; j < mvz; j++) {
      wald[(2*i+1)*nWald + 4 + j] = REAL(est)[j] + sol_w*REAL(hFunc -> vz)[j];
    }    
    /* keep track of status */
    if ((i+1) % 200 == 0 && i != 0 && i != mz - 1) {
      Rprintf("worker pid=%d generated %d rays\n", (int) getpid(), i+1);
    }
  }
  Rprintf("Worker pid=%d is done generating %d rays!\n", (int) getpid(), i);
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
  Free(out);
  Free(wald);
  UNPROTECT(16);
  return output;
}
