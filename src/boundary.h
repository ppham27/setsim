
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

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


#endif /* BOUNDARY_H_ */
