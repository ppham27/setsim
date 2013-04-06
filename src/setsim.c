#include "setsim.h"

/* get the list element named str, or return NULL */
int getListIndex(SEXP list, const char *str) {
  SEXP names = getAttrib(list, R_NamesSymbol);  
  for (int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      return i;
    }
  }
  Rprintf("Error: %s not found in list\n", str);
  return -1;
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
