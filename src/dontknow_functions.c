#define USE_FC_LEN_T

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
// #include <omp.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rconfig.h>

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#ifndef FCONE
# define FCONE
#endif

/* donknow family log-likelihood. */
/* clamp helper */
static inline double clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

/* Bivariate normal CDF Φ2(h,k;rho) over (-Inf,h]×(-Inf,k] */
static double bvncdf(double h, double k, double rho)
{
  /* Guard rho strictly inside (-1,1) */
  if (rho >= 1.0)  rho =  0.999999999;
  if (rho <= -1.0) rho = -0.999999999;

  /* Quick exits for infinities */
  if (!R_FINITE(h)) {
    if (h < 0) return 0.0;
    return pnorm(k, 0.0, 1.0, 1, 0);
  }
  if (!R_FINITE(k)) {
    if (k < 0) return 0.0;
    return pnorm(h, 0.0, 1.0, 1, 0);
  }

  /* Nearly independent? Use product */
  if (fabs(rho) < 1e-12) {
    return pnorm(h, 0.0, 1.0, 1, 0) * pnorm(k, 0.0, 1.0, 1, 0);
  }

  /* 5-point Gauss–Legendre nodes/weights on [-1,1] */
  static const double w[5] = {
    0.017614007139152118, 0.040601429800386941, 0.062672048334109064,
    0.083276741576704749, 0.10193011981724044
  };
  static const double x[5] = {
    0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
    0.8391169718222188, 0.7463319064601508
  };

  const double h2 = h*h, k2 = k*k;

  /* Two-regime algorithm (Genz/Drezner): */
  double result;
  if (fabs(rho) < 0.925) {
    const double asr = asin(rho);
    double sum = 0.0;
    for (int i = 0; i < 5; ++i) {
      for (int s = -1; s <= 1; s += 2) {
        double sn = sin(asr * (1.0 + s * x[i]) * 0.5);
        double ss = 1.0 - sn*sn;
        double ex = exp( (2.0*sn*h*k - h2 - k2) / (2.0*ss) ) / sqrt(ss);
        sum += w[i] * ex;
      }
    }
    result = pnorm(h,0,1,1,0) * pnorm(k,0,1,1,0) + (asr * sum) / (2.0 * M_PI);
  } else {
    /* Near |rho| ~ 1: alternative stable form */
    const double rr = 1.0 - rho*rho;
    const double s  = (rho < 0.0) ? -1.0 : 1.0;

    /* Leading term */
    double a = (k - s*h) / sqrt(rr);
    double term1 = 0.0;
    if (s*h <= k) term1 = pnorm(a, 0.0, 1.0, 1, 0);

    /* Correction integral */
    double sum = 0.0;
    double hk  = h*k;
    for (int i = 0; i < 5; ++i) {
      for (int sgn = -1; sgn <= 1; sgn += 2) {
        double r = (1.0 + sgn * x[i]) * 0.5;    /* map to [0,1] */
        double t2 = 1.0 - r * rr;               /* t^2 */
        double ex = exp( -(h2 - 2.0*rho*hk + k2) / (2.0 * (1.0 - t2)) ) / (1.0 - t2);
        sum += w[i] * ex;
      }
    }
    result = term1 * exp(-0.5*(h2+k2)) / (2.0*M_PI) + s * sum / (2.0*M_PI);
    result = clamp(result, 0.0, 1.0);
  }
  return clamp(result, 0.0, 1.0);
}

/* .Call entry:
   logLik_dontknow(eta1, eta2, rho, alpha, y, logOut)
   - eta1, eta2: REAL vectors length n
   - rho: REAL vector length 1 or n (correlation scale)
   - alpha: REAL matrix n×4 (columns: alpha1,alpha2,alpha3,alpha4)
   - y: REAL/INT matrix n×2 (y1,y2)
   - logOut: logical
*/
SEXP logLik_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha, SEXP Y, SEXP LogOut)
{
  if (!isReal(Eta1) || !isReal(Eta2)) error("eta1/eta2 must be real.");
  if (!isReal(Rho)) error("rho must be real.");
  if (!isMatrix(Alpha) || !isReal(Alpha)) error("alpha must be a real matrix n x 4.");
  if (!isMatrix(Y)) error("y must be a matrix n x 2.");

  const R_xlen_t n = XLENGTH(Eta1);
  if (XLENGTH(Eta2) != n) error("eta1 and eta2 must have same length.");

  SEXP dimA = getAttrib(Alpha, R_DimSymbol);
  if (TYPEOF(dimA) != INTSXP || INTEGER(dimA)[0] != n || INTEGER(dimA)[1] != 4)
    error("alpha must have dimensions n x 4.");

  SEXP dimY = getAttrib(Y, R_DimSymbol);
  if (TYPEOF(dimY) != INTSXP || INTEGER(dimY)[0] != n || INTEGER(dimY)[1] != 2)
    error("y must have dimensions n x 2.");

  const double *eta1 = REAL(Eta1);
  const double *eta2 = REAL(Eta2);
  const double *rhoV = REAL(Rho);
  const double *A    = REAL(Alpha);
  const double *Yp   = REAL(Y);      /* works if Y is integer or real */

  const int rho_scalar = (XLENGTH(Rho) == 1);
  int log_out = asLogical(LogOut);
  if (log_out == NA_LOGICAL) log_out = 1;

  SEXP ans = PROTECT(allocVector(REALSXP, n));
  double *ll = REAL(ans);

  for (R_xlen_t i = 0; i < n; ++i) {
    /* row-major helpers */
    const double a1 = A[i + 0*n];  /* alpha1(i) */
    const double a2 = A[i + 1*n];
    const double a3 = A[i + 2*n];
    const double a4 = A[i + 3*n];

    /* y1,y2 */
    const int y1 = (int) Yp[i + 0*n];
    const int y2 = (int) Yp[i + 1*n];

    /* rho for obs i */
    double r = rho_scalar ? rhoV[0] : rhoV[i];
    r = clamp(r, -0.999999, 0.999999);

    if (y1 == 1) {
      /* dk: P(Y1=1) = 1 - Phi(alpha1 - eta1) */
      const double z = a1 - eta1[i];
      double p = pnorm(z, 0.0, 1.0, 0, 0);  /* upper tail = 1 - Phi(z) */
      /* numeric guard */
      if (p < 0.0) p = 0.0;
      if (log_out) ll[i] = (p > 0.0) ? log(p) : R_NegInf;
      else         ll[i] = p;
    } else {
      /* y1 == 0: rectangle */
      const double Aup = a1 - eta1[i];   /* first-dim upper bound */

      /* build second-dim full cuts: (-Inf, a2, a3, a4, +Inf) */
      /* choose bounds by y2 in {0,1,2,3} */
      double Bl, Bu;
      if (y2 == 0) {        Bl = R_NegInf;    Bu = a2 - eta2[i]; }
      else if (y2 == 1) {   Bl = a2 - eta2[i]; Bu = a3 - eta2[i]; }
      else if (y2 == 2) {   Bl = a3 - eta2[i]; Bu = a4 - eta2[i]; }
      else { /* y2 == 3 */  Bl = a4 - eta2[i]; Bu = R_PosInf; }

      double p_up = bvncdf(Aup, Bu, r);
      double p_lo = bvncdf(Aup, Bl, r);
      double p = p_up - p_lo;

      /* guards */
      if (p < 0.0) p = 0.0;
      if (log_out) ll[i] = (p > 0.0) ? log(p) : R_NegInf;
      else         ll[i] = p;
    }
  }

  UNPROTECT(1);
  return ans;
}

