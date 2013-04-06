#include "setsim.h"
SEXP runUnitTests(SEXP h, SEXP hEnv, SEXP input) {
  SEXP successful;
  PROTECT(successful = allocVector(LGLSXP, 1));
  LOGICAL(successful)[0] = 1;

  /* test VZ matrix */
  int mv = 3;
  int nv = 3;
  double v[9] = {0.15306913, -0.05906771, -0.07683781,
                 -0.059067706, 0.154188794, 0.003284563,
                 -0.076837810, 0.003284563, 0.158885972};
  int mz = 2;
  int nz = 3;
  double z[6] = {-1.123198813, -1.919282365,
                 -0.148169045, 1.276729453,
                 -0.357214516, -0.509841023};
  double vz_correct[6] = {-0.13572747, 0.04232548, 0.02906109,
                          -0.33002129, 0.30855038, 0.07066037};
  /* calculate VZ matrix */
  double *vz_test = calculateVZMatrix(mv, nv, v,
                                      mz, nz, z);
  Rprintf("VZ matrix evaluation test:\n");
  for (int i = 0; i < mv*mz; i++) {
    if ( fabs(vz_correct[i] - vz_test[i]) > 0.000001 ) {
      Rprintf("%f does not match %f\n",vz_test[i],vz_correct[i]);
      LOGICAL(successful)[0] = 0;
    }
  }
  if (LOGICAL(successful)[0] == 1) {
    Rprintf("VZ matrix evaluation was successful.\n");
  }

  /* test hFunction evaluation */
  hFunction *poisson;
  poisson = (hFunction *) R_alloc(1, sizeof(hFunction));
  poisson -> h = h;
  poisson -> env = hEnv;
  poisson -> epsilon = Rf_ScalarReal(-60.45041);
  poisson -> input = input;
  SEXP vz_poisson;
  vz_poisson = allocVector(REALSXP, 3);
  REAL(vz_poisson)[0] = 0.01067525;
  REAL(vz_poisson)[1] = -0.01109185;
  REAL(vz_poisson)[2] = -0.01262752;
  poisson -> vz = vz_poisson;
  poisson -> cstar = Rf_ScalarReal(2.365974);
  poisson -> fcall = lang5(poisson -> h,
                           poisson -> epsilon,
                           poisson -> input,
                           poisson -> vz,
                           poisson -> cstar);

  double h_correct[2] = {-22.07824, -14.152693};
  double *h_test;
  h_test = (double *) R_alloc(2, sizeof(double));
  h_test[0] = REAL(eval(poisson -> fcall,
                        poisson -> env))[0];
  /* change inputs */
  REAL(poisson -> epsilon)[0] = -42.97172;
  REAL(poisson -> vz)[0] = 0.015725877;
  REAL(poisson -> vz)[1] = -0.019692698;
  REAL(poisson -> vz)[2] = 0.001130277;
  REAL(poisson -> cstar)[0] = 2.365974;
  h_test[1] = REAL(eval(poisson -> fcall,
                        poisson -> env))[0];
  Rprintf("H function evaluation test:\n");
  for (int i = 0; i < 2; i++) {
    if ( fabs(h_correct[i] - h_test[i]) > 0.000001 ) {
      Rprintf("%f does not match %f\n",h_test[i],h_correct[i]);
      LOGICAL(successful)[0] = 0;
    }
  }
  if (LOGICAL(successful)[0] == 1) {
    Rprintf("H function evaluation was successful.\n");
  }

  /* test nroots */
  int nhValues = 51;
  double hValues[51] = {-21.5287242753573,-19.5890883182986,-17.7497401154925,
                        -16.0085795506889,-14.3635517673602,-12.8126461678268,
                        -11.3538954350094,-9.98537457629332,-8.70519998899077,
                        -7.51152854690595,-6.40255670751814,-5.37651963930715,
                        -4.43169036875885,-3.56637894659884,-2.77893163281173,
                        -2.06773010001474,-1.43119065476581,-0.86776347639106,
                        -0.375931872935196,0.0457884461644524,0.39885008105678
                        ,0.684674636488596,0.904653401212401,1.06014801214722,
                        1.1524911036429,1.18298694218767,1.15291204689199,
                        1.06351579607453,0.916021020266997,0.711624581948074,
                        0.451497942310952,0.136787715358496,-0.23138379038373,
                        -0.651918042254251,-1.12373976566648,-1.64579643787441,
                        -2.21705779302947,-2.83651533827144,-3.50318188060261,
                        -4.21609106429993,-4.97429691862426,-5.776873415595,
                        -6.62291403759822,-7.51153135460767,-8.44185661079909,
                        -9.41303932034439,-10.4242468721783,-11.4746641435332,
                        -12.5634931220439,-13.6899525362276,-14.85327749415};
  int nepsilon = 51;
  double epsilon[51] = {-45.6629511583333,-43.836433112,-42.0099150656666,
                        -40.1833970193333,-38.356878973,-36.5303609266666,
                        -34.7038428803333,-32.877324834,-31.0508067876666,
                        -29.2242887413333,-27.397770695,-25.5712526486667,
                        -23.7447346023333,-21.918216556,-20.0916985096667,
                        -18.2651804633333,-16.438662417,-14.6121443706667,
                        -12.7856263243333,-10.959108278,-9.13259023166666,
                        -7.30607218533333,-5.479554139,-3.65303609266667,
                        -1.82651804633333,0,1.82651804633333,
                        3.65303609266666,5.47955413899999,7.30607218533333,
                        9.13259023166665,10.959108278,12.7856263243333,
                        14.6121443706667,16.438662417,18.2651804633333,
                        20.0916985096667,21.918216556,23.7447346023333,
                        25.5712526486666,27.397770695,29.2242887413333,
                        31.0508067876666,32.877324834,34.7038428803333,
                        36.5303609266666,38.356878973,40.1833970193333,
                        42.0099150656666,43.836433112,45.6629511583333};
  int nnroots = 2;
  int index[2] = {19, 32};
  /* Browse[2]> sol */
  /*            x.lower   x.upper    y.lower     y.upper */
  /* Interval -12.78563 -10.95911 -0.3759319  0.04578845 */
  /* Interval  10.95911  12.78563  0.1367877 -0.23138379 */
  double sol[8] = {-12.78563,-10.95911,-0.3759319,0.04578845,
                   10.95911,12.78563,0.1367877,-0.23138379};
  int *position;
  double *range;
  position = (int *) R_alloc(nhValues, sizeof(int));
  range = (double *) R_alloc(nhValues*4, sizeof(double));
  int n = nroots(nhValues, hValues,
                 nepsilon, epsilon,
                 position,
                 range);
  Rprintf("nroots evaluation test:\n");
  if ( n != nnroots ) {
    Rprintf("%d does not match %d\n", n, nnroots);
    LOGICAL(successful)[0] = 0;
  }
  for (int i = 0; i < n; i++) {
    if (index[i] != position[i]) {
      Rprintf("%d does not match %d\n", position[i], index[i]);
      LOGICAL(successful)[0] = 0;
    }
  }
  for (int i = 0; i < n*4; i++) {
    if ( fabs(range[i] - sol[i]) > 0.00001 ) {
      Rprintf("%f does not match %f\n", range[i], sol[i]);
      LOGICAL(successful)[0] = 0;
    }
  }
  if (LOGICAL(successful)[0] == 1) {
    Rprintf("nroots evaluation was successful.\n");
  }


  REAL(poisson -> vz)[0] = 0.013510572;
  REAL(poisson -> vz)[1] = -0.012675022;
  REAL(poisson -> vz)[2] = 0.004541779;
  Rprintf("bisect evaluation test:\n");
  double testRoot = bisect(poisson, -11.53861, -9.890235);
  double realRoot = -10.604734;
  if ( fabs(realRoot - testRoot) > 0.00001 ) {
    Rprintf("%f does not match %f\n", testRoot, realRoot);
    LOGICAL(successful)[0] = 0;
  } else {
    Rprintf("bisect evaluation was successful.\n");
  }
  

  UNPROTECT(1);
  return successful;
}


