// src/logLik_dontknow.c
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/*
  Bivariate normal CDF Phi2(h, k; rho) using Plackett's integral:
  Phi2(h,k;rho) = Phi(h) * Phi(k) + ∫_{0}^{rho} f(t) dt,
  where f(t) = [1/(2π sqrt(1 - t^2))] * exp( -(h^2 - 2 h k t + k^2) / (2(1 - t^2)) ).
  We evaluate the integral with 20-point Gauss–Legendre on [min(0,rho), max(0,rho)].
  For rho=0 we return Phi(h)Phi(k). Accurate for |rho|<1 (typical statistical use).
*/

static double Phi(double x) {
  // lower_tail=1 (TRUE), log_p=0 (FALSE)
  return Rf_pnorm5(x, 0.0, 1.0, 1, 0);
}

static double phi2_integrand(double h, double k, double t) {
  double one_minus_t2 = 1.0 - t * t;
  double denom = 2.0 * one_minus_t2;
  double quad = (h*h - 2.0*h*k*t + k*k) / denom;
  double root = sqrt(one_minus_t2);
  // 1/(2π sqrt(1 - t^2)) * exp( -quad )
  return exp(-quad) / (2.0 * M_PI * root);
}

// 20-point Gauss–Legendre nodes and weights on [-1, 1]
static const double GL20_x[20] = {
  -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188, -0.7463319064601508,
  -0.6360536807265150, -0.5108670019508271, -0.3737060887154196, -0.2277858511416451, -0.07652652113349733,
   0.07652652113349733,  0.2277858511416451,  0.3737060887154196,  0.5108670019508271,  0.6360536807265150,
   0.7463319064601508,  0.8391169718222188,  0.9122344282513259,  0.9639719272779138,  0.9931285991850949
};

static const double GL20_w[20] = {
  0.01761400713915212, 0.04060142980038694, 0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
  0.1181945319615184,  0.1316886384491766,  0.1420961093183821,  0.1491729864726037,  0.1527533871307259,
  0.1527533871307259,  0.1491729864726037,  0.1420961093183821,  0.1316886384491766,  0.1181945319615184,
  0.1019301198172404,  0.08327674157670475, 0.06267204833410906, 0.04060142980038694, 0.01761400713915212
};

static double Phi2(double h, double k, double rho)
{
  if (rho == 0.0) return Phi(h) * Phi(k);

  // Integrate from a = min(0, rho) to b = max(0, rho)
  double a = (rho < 0.0) ? rho : 0.0;
  double b = (rho < 0.0) ? 0.0 : rho;

  // Affine map from [-1,1] to [a,b]: t = ((b-a)/2)*x + (a+b)/2
  double mid = 0.5 * (a + b);
  double half = 0.5 * (b - a);

  double sum = 0.0;
  for (int i = 0; i < 20; ++i) {
    double t = half * GL20_x[i] + mid;
    sum += GL20_w[i] * phi2_integrand(h, k, t);
  }

  double integral = half * sum;
  return Phi(h) * Phi(k) + integral;
}

/*
  .Call interface:
    logLik_dontknow(eta1, eta2, rho, alpha, y, log)
  Arguments:
    - eta1: numeric vector (length n)
    - eta2: numeric vector (length n)
    - rho:  numeric scalar or vector (recycled to length n)
    - alpha: numeric matrix n x k (cols: alpha1,...,alphak), with alpha[,-1] already converted to cumulative cuts
    - y: integer/numeric matrix n x 2 (cols: y1 ∈ {0,1}, y2 ∈ {0,1,...,k-1})
    - log: logical scalar
  Returns:
    - numeric vector length n with log-likelihood contributions if log=TRUE, else probabilities.
*/

SEXP logLik_dontknow(SEXP ETA1, SEXP ETA2, SEXP RHO, SEXP ALPHA, SEXP Y, SEXP LOGP)
{
  if (!isReal(ETA1) || !isReal(ETA2)) error("eta1 and eta2 must be numeric (double).");
  if (!isReal(RHO)) error("rho must be numeric (double).");
  if (!isReal(ALPHA) || !isMatrix(ALPHA)) error("alpha must be a numeric matrix.");
  if ((!isInteger(Y) && !isReal(Y)) || !isMatrix(Y)) error("y must be an integer/numeric matrix with two columns.");
  if (!isLogical(LOGP) || LENGTH(LOGP) != 1) error("log must be a single logical.");

  R_xlen_t n = XLENGTH(ETA1);
  if (XLENGTH(ETA2) != n) error("eta1 and eta2 must have the same length.");

  SEXP dimA = getAttrib(ALPHA, R_DimSymbol);
  if (isNull(dimA) || LENGTH(dimA) != 2) error("alpha must have dimensions.");
  R_xlen_t nA = (R_xlen_t) INTEGER(dimA)[0];
  int k = INTEGER(dimA)[1];
  if (nA != n) error("nrow(alpha) must match length of eta1.");

  SEXP dimY = getAttrib(Y, R_DimSymbol);
  if (isNull(dimY) || LENGTH(dimY) != 2) error("y must have dimensions.");
  R_xlen_t nY = (R_xlen_t) INTEGER(dimY)[0];
  int cy = INTEGER(dimY)[1];
  if (nY != n || cy != 2) error("y must be an n x 2 matrix.");

  const double *eta1 = REAL(ETA1);
  const double *eta2 = REAL(ETA2);
  const double *rho  = REAL(RHO);
  const double *alpha = REAL(ALPHA);

  // y can be integer or real — access via REAL coercion when needed
  int y_is_int = isInteger(Y);
  const int   *yI = y_is_int ? INTEGER(Y) : NULL;
  const double*yR = y_is_int ? NULL : REAL(Y);

  int logp = LOGICAL(LOGP)[0];

  SEXP ans = PROTECT(allocVector(REALSXP, n));
  double *out = REAL(ans);

  for (R_xlen_t i = 0; i < n; ++i) {
    // Observed (y1, y2)
    int y1 = y_is_int ? (int) yI[i] : (int) floor(yR[i] + 0.5);                 // column 1
    int y2 = y_is_int ? (int) yI[i + n] : (int) floor(yR[i + n] + 0.5);         // column 2

    // Correlation (recycle if needed)
    double r = REAL(RHO)[ (XLENGTH(RHO) == 1) ? 0 : i ];
    if (r <= -0.999999) r = -0.999999;
    if (r >=  0.999999) r =  0.999999;

    // alpha row i: alpha[i, 0..k-1]
    const double *ai = alpha + i + 0 * n; // start of col 0, row i
    // helper to index alpha[i, j] = alpha[i + j*n]
    #define A(j) (*(ai + (size_t)(j) * (size_t)n))

    double A1 = A(0);                       // alpha1
    double H  = A1 - eta1[i];               // A = alpha1 - eta1

    if (y1 == 1) {
      // P(Y1=1) = 1 - Phi(alpha1 - eta1)
      double val = Rf_pnorm5(H, 0.0, 1.0, 0, logp); // lower_tail=FALSE
      out[i] = val;
    } else {
      // y1 == 0: rectangle probability difference on (Z1, Z2)
      // Build implicit alpha2 with -Inf and +Inf around alpha[,-1]
      // alpha2 indices: 0 -> -Inf, 1..(k-1) -> alpha[,1..k-1], k -> +Inf
      int c_lo = y2 + 1;     // index in alpha2 for lower bound
      int c_up = y2 + 2;     // index in alpha2 for upper bound

      double B1, B2;
      if (c_lo == 0) {
        B1 = R_NegInf;
      } else {
        // alpha[, c_lo] corresponds to alpha column (c_lo) in 1-based,
        // i.e., A(c_lo-1) in 0-based since alpha[,-1] shifted by one.
        B1 = A(c_lo - 1);
      }
      if (c_up >= k + 1) {
        B2 = R_PosInf;
      } else {
        B2 = A(c_up - 1);
      }

      // Shift by eta2
      B1 -= eta2[i];
      B2 -= eta2[i];

      // Evaluate P = Phi2(H, B2; r) - Phi2(H, B1; r)
      double pup = Phi2(H, B2, r);
      double plo = Phi2(H, B1, r);
      double p = pup - plo;

      if (logp) {
        // Guard tiny numerical negatives
        out[i] = (p > 0.0) ? log(p) : R_NegInf;
      } else {
        out[i] = (p > 0.0) ? p : 0.0;
      }
    }
    #undef A
  }

  UNPROTECT(1);
  return ans;
}

