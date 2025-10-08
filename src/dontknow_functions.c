#define USE_FC_LEN_T

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* ---- helpers ---- */
static inline double clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

/* Deterministic bivariate normal CDF Φ2(h,k;rho) over (-Inf,h]×(-Inf,k] */
static double bvncdf(double h, double k, double rho)
{
  if(rho >= 1.0)  rho =  0.999999999;
  if(rho <= -1.0) rho = -0.999999999;

  if(!R_FINITE(h)) return (h < 0.0) ? 0.0 : pnorm(k, 0.0, 1.0, 1, 0);
  if(!R_FINITE(k)) return (k < 0.0) ? 0.0 : pnorm(h, 0.0, 1.0, 1, 0);

  if(fabs(rho) < 1e-12) {
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
  double result;

  if(fabs(rho) < 0.925) {
    const double asr = asin(rho);
    double sum = 0.0;
    for(int i = 0; i < 5; ++i) {
      for(int s = -1; s <= 1; s += 2) {
        const double sn = sin(asr * (1.0 + s * x[i]) * 0.5);
        const double ss = 1.0 - sn*sn;
        const double ex = exp((2.0*sn*h*k - h2 - k2) / (2.0*ss)) / sqrt(ss);
        sum += w[i] * ex;
      }
    }
    result = pnorm(h,0,1,1,0) * pnorm(k,0,1,1,0) + (asr * sum) / (2.0 * M_PI);
  } else {
    const double rr = 1.0 - rho*rho;
    const double s  = (rho < 0.0) ? -1.0 : 1.0;

    const double a = (k - s*h) / sqrt(rr);
    const double lead = (s*h <= k) ? pnorm(a, 0.0, 1.0, 1, 0) : 0.0;

    double sum = 0.0, hk = h*k;
    for(int i = 0; i < 5; ++i) {
      for(int sg = -1; sg <= 1; sg += 2) {
        const double r = (1.0 + sg * x[i]) * 0.5;
        const double denom = 1.0 - (1.0 - rho*rho) * r;
        const double ex = exp(-(h2 - 2.0*rho*hk + k2) / (2.0 * denom)) / denom;
        sum += w[i] * ex;
      }
    }
    result = lead * exp(-0.5*(h2+k2)) / (2.0*M_PI) + s * sum / (2.0*M_PI);
    result = clamp(result, 0.0, 1.0);
  }
  return clamp(result, 0.0, 1.0);
}

/* .Call interface:
   logLik_dontknow(eta1, eta2, rho, alpha, y, logOut)

   - eta1, eta2: REAL vectors length n
   - rho: REAL vector length 1 or n (correlation scale)
   - alpha: REAL matrix n×m (m>=2), cols: [alpha1 | alpha2..alpha_m], ordered row-wise for cols>=2
   - y: INT/REAL matrix n×2 (y1∈{0,1}, y2∈{0,..,m-2})
   - logOut: logical
*/
SEXP logLik_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha, SEXP Y, SEXP LogOut)
{
  if(!isReal(Eta1) || !isReal(Eta2)) error("eta1/eta2 must be real.");
  if(!isReal(Rho)) error("rho must be real.");
  if(!isMatrix(Alpha) || !isReal(Alpha)) error("alpha must be real matrix n x m (m>=2).");
  if(!isMatrix(Y)) error("y must be matrix n x 2.");

  const R_xlen_t n = XLENGTH(Eta1);
  if(XLENGTH(Eta2) != n) error("eta1 and eta2 must have same length.");

  SEXP dimA = getAttrib(Alpha, R_DimSymbol);
  if(TYPEOF(dimA) != INTSXP) error("alpha must have integer dims.");
  const int nA = INTEGER(dimA)[0];
  const int mA = INTEGER(dimA)[1];
  if(nA != (int)n || mA < 2) error("alpha must have dimensions n x m with m>=2.");

  SEXP dimY = getAttrib(Y, R_DimSymbol);
  if(TYPEOF(dimY) != INTSXP) error("y must have integer dims.");
  const int nY = INTEGER(dimY)[0];
  const int mY = INTEGER(dimY)[1];
  if(nY != (int)n || mY != 2) error("y must have dimensions n x 2.");

  const double *eta1 = REAL(Eta1);
  const double *eta2 = REAL(Eta2);
  const double *rhoV = REAL(Rho);
  const double *A    = REAL(Alpha);   /* col-major: A[i + j*n] */
  const double *Yp   = REAL(Y);

  const int rho_scalar = (XLENGTH(Rho) == 1);
  int log_out = asLogical(LogOut);
  if(log_out == NA_LOGICAL) log_out = 1;

  /* number of ordinal levels in Y2 */
  const int K = mA - 1;         /* categories for Y2 are 0..K-1 */

  SEXP out = PROTECT(allocVector(REALSXP, n));
  double *ll = REAL(out);

  for(R_xlen_t i = 0; i < n; ++i) {
    const double a1 = A[i + 0*n];      /* binary cut alpha1(i) */

    const int y1 = (int) Yp[i + 0*n];
    int y2       = (int) Yp[i + 1*n];

    /* basic checks on y1/y2 */
    if(y1 != 0 && y1 != 1) {
      error("y1 must be 0 or 1 at row %lld.", (long long)(i+1));
    }
    if(y2 < 0 || y2 >= K) {
      error("y2 out of range [0,%d] at row %lld.", K-1, (long long)(i+1));
    }

    double r = rho_scalar ? rhoV[0] : rhoV[i];
    r = clamp(r, -0.999999, 0.999999);

    if(y1 == 1) {
      /* P(Y1=1) = 1 - Phi(alpha1 - eta1) */
      const double z = a1 - eta1[i];
      double p = pnorm(z, 0.0, 1.0, 0, 0);
      if(p < 0.0) p = 0.0;
      ll[i] = log_out ? ((p > 0.0) ? log(p) : R_NegInf) : p;
    } else {
      /* y1==0: rectangle on (Z1, Z2) */
      const double Aup = a1 - eta1[i];

      /* Build lower/upper bounds for Z2 based on y2 category */
      double Bl, Bu;
      if(y2 == 0) {
        Bl = R_NegInf;
        Bu = A[i + 1*n] - eta2[i];                 /* alpha2 */
      } else if(y2 == K-1) {
        Bl = A[i + (K)*n] - eta2[i];               /* alpha_{K+1} column index = K */
        Bu = R_PosInf;
      } else {
        Bl = A[i + (y2+1)*n] - eta2[i];            /* alpha_{y2+1} */
        Bu = A[i + (y2+2)*n] - eta2[i];            /* alpha_{y2+2} */
      }

      double p = bvncdf(Aup, Bu, r) - bvncdf(Aup, Bl, r);
      if(p < 0.0) p = 0.0;                        /* guard tiny negatives */
      ll[i] = log_out ? ((p > 0.0) ? log(p) : R_NegInf) : p;
    }
  }

  UNPROTECT(1);
  return out;
}

