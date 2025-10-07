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

/* bivariate normal CDF over (-inf,h] x (-inf,k] with correlation rho */
static double bvncdf(double h, double k, double rho)
{
  /* handle extreme correlations safely */
  if(rho >= 1.0)  rho =  0.999999999;
  if(rho <= -1.0) rho = -0.999999999;

  /* trivial cases */
  if(!R_FINITE(h)) {
    if(h < 0) return 0.0;
    /* h = +Inf */
    return pnorm(k, 0.0, 1.0, 1, 0);
  }
  if(!R_FINITE(k)) {
    if(k < 0) return 0.0;
    /* k = +Inf */
    return pnorm(h, 0.0, 1.0, 1, 0);
  }

  /* if nearly independent, fall back to product */
  if(fabs(rho) < 1e-12) {
    return pnorm(h, 0.0, 1.0, 1, 0) * pnorm(k, 0.0, 1.0, 1, 0);
  }

  /* Constants / nodes-weights for Gauss-Legendre */
  static const double w[5] = {
      0.017614007139152118, 0.040601429800386941, 0.062672048334109064,
      0.083276741576704749, 0.10193011981724044
  };
  static const double x[5] = {
      0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
      0.8391169718222188, 0.7463319064601508
  };

  const double ph = pnorm(h, 0.0, 1.0, 1, 0);
  const double pk = pnorm(k, 0.0, 1.0, 1, 0);

  const double rr = 1.0 - rho*rho;
  const double h2 = h*h, k2 = k*k;

  double result;

  /* Two-regime algorithm:
     - |rho| < 0.925: Genz/Drezner formula using asr = asin(rho)
     - otherwise: alternative stable form near |rho| ~ 1
  */
  if(fabs(rho) < 0.925) {
    double sum = 0.0;
    const double asr = asin(rho);

    for(int i = 0; i < 5; ++i) {
      for(int sgn = -1; sgn <= 1; sgn += 2) {
        double sn = sin(asr * (1.0 + sgn * x[i]) * 0.5);
        double ss = 1.0 - sn*sn;
        double ex = exp( (2.0*sn*h*k - h2 - k2) / (2.0*ss) ) / sqrt(ss);
        sum += w[i] * ex;
      }
    }
    result = ph*pk + (asr * sum) / (2.0 * M_PI);
  } else {
    /* Near |rho| ~ 1, use alternative expression */
    double s = (rho < 0.0) ? -1.0 : 1.0;
    double a = (k - s*h) / sqrt(rr);
    double b = (k + s*h) / sqrt(rr);

    /* Tail-safe handling */
    double term1 = 0.0;
    if(s*h <= k) {
      term1 = pnorm(a, 0.0, 1.0, 1, 0);
    }

    /* correction integral using Gauss-Legendre */
    double sum = 0.0;
    double hk = h*k;
    for(int i = 0; i < 5; ++i) {
      for(int sgn = -1; sgn <= 1; sgn += 2) {
        double r = (1.0 + sgn * x[i]) / 2.0;
        double t = sqrt(1.0 - r * rr);
        double ex = exp( -(h2 - 2.0*rho*hk + k2) / (2.0*(1.0 - t*t)) ) / (1.0 - t*t);
        sum += w[i] * ex;
      }
    }
    result = term1 * exp( -0.5 * (h2 + k2) ) / (2.0 * M_PI) + s * sum / (2.0 * M_PI);
    /* transform back to CDF scale */
    result = clamp(result, 0.0, 1.0); /* numerical guard */
    /* This branch is a protective approximation; refine with product baseline */
    if(rho > 0) result = fmax(result, ph*pk); /* keep monotone-ish */
  }

  /* final guard */
  result = clamp(result, 0.0, 1.0);
  return result;
}

SEXP logLik_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha1, SEXP Alpha2, SEXP Y, SEXP LogOut)
{
  if(!isReal(Eta1) || !isReal(Eta2) || !isReal(Rho) ||
      !isReal(Alpha1) || !isReal(Alpha2) || !isMatrix(Y))
    error("Invalid input types.");

  R_xlen_t n = XLENGTH(Eta1);
  if(XLENGTH(Eta2) != n) error("eta1 and eta2 must have same length.");

  SEXP dimY = getAttrib(Y, R_DimSymbol);
  if(TYPEOF(dimY) != INTSXP || INTEGER(dimY)[1] != 2)
    error("y must be an n x 2 matrix.");

  const double *eta1 = REAL(Eta1);
  const double *eta2 = REAL(Eta2);
  const double *rhoV = REAL(Rho);
  const double *a1   = REAL(Alpha1);
  const double *a2   = REAL(Alpha2);
  const double *y    = REAL(Y); /* allow integer-coerced */

  int log_out = asLogical(LogOut);
  if(log_out == NA_LOGICAL) log_out = 1;

  int rho_scalar = (XLENGTH(Rho) == 1);

  SEXP ans = PROTECT(allocVector(REALSXP, n));
  double *ll = REAL(ans);

  for(R_xlen_t i = 0; i < n; ++i) {
    /* y1, y2 are coded {0,1,...} as in your R code */
    int y1 = (int) y[i + 0*n]; /* first column */
    int y2 = (int) y[i + 1*n]; /* second column */
    double r  = rho_scalar ? rhoV[0] : rhoV[i];
    r = clamp(r, -0.999999, 0.999999);

    /* thresholds:
       alpha1[y1+1]   and alpha1[y1+2] exist because alpha1 = (-Inf, a, +Inf)
       alpha2[y2+1]   and alpha2[y2+2] for the ordinal part
    */
    double u1 = a1[y1 + 1]; /* upper on dim 1 */
    double l1 = (y1 == 0) ? R_NegInf : a1[y1]; /* lower on dim 1; here we only need (-inf, u1] */

    if(y1 == 1) {
      /* dk: marginal on dim 1 only: P(a1[y1] < Z1 <= a1[y1+1]) = Phi(u1-eta1) - Phi(l1-eta1) */
      double p1 = pnorm(u1 - eta1[i], 0.0, 1.0, 1, 0) - pnorm(a1[y1] - eta1[i], 0.0, 1.0, 1, 0);
      p1 = fmax(p1, 0.0); /* guard */
      ll[i] = log_out ? ( (p1 > 0.0) ? log(p1) : R_NegInf ) : p1;
    } else {
      /* non-dk: rectangle difference on second margin */
      double A = u1 - eta1[i];

      double B2u = a2[y2 + 2] - eta2[i]; /* upper bound for Y2 */
      double B1u = a2[y2 + 1] - eta2[i]; /* lower bound for Y2 */

      /* P(Z1<=A, Z2<=B2u) - P(Z1<=A, Z2<=B1u) */
      double pA_B2 = bvncdf(A, B2u, r);
      double pA_B1 = bvncdf(A, B1u, r);

      double p = pA_B2 - pA_B1;
      if(p < 0.0) p = 0.0; /* guard small negatives */
      ll[i] = log_out ? ( (p > 0.0) ? log(p) : R_NegInf ) : p;
    }
  }

  UNPROTECT(1);
  return ans;
}

