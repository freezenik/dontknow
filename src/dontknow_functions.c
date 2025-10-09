#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* If you have a miwa.h, include it; otherwise forward-declare R_miwa */
SEXP R_miwa(SEXP steps, SEXP corr, SEXP upper, SEXP lower, SEXP infin);

/* helper: single BVN cdf Phi_2( (-Inf, -Inf) -> (h, k) ; rho ) via R_miwa */
static inline double bvn_cdf_miwa(double h, double k, double rho, int steps)
{
    /* Lower bounds (-Inf, -Inf), both upper-bounded (code 0) */
    SEXP Steps = PROTECT(ScalarInteger(steps));

    /* 2x2 correlation (column-major doesn’t matter; R_miwa reads symmetric) */
    SEXP Corr  = PROTECT(allocVector(REALSXP, 4));
    double *c  = REAL(Corr);
    c[0] = 1.0;  c[1] = rho;
    c[2] = rho;  c[3] = 1.0;

    SEXP Upper = PROTECT(allocVector(REALSXP, 2));
    double *u  = REAL(Upper);
    u[0] = h;   u[1] = k;

    SEXP Lower = PROTECT(allocVector(REALSXP, 2));
    double *lo = REAL(Lower);
    lo[0] = R_NegInf; lo[1] = R_NegInf;

    SEXP Infin = PROTECT(allocVector(INTSXP, 2));
    int *inf   = INTEGER(Infin);
    inf[0] = 0; inf[1] = 0; /* (-Inf, upper] for both dims */

    SEXP res = PROTECT(R_miwa(Steps, Corr, Upper, Lower, Infin));
    double ans = REAL(res)[0];

    UNPROTECT(6);
    return ans;
}

/* .Call entry:
     logLik_dontknow(eta1, eta2, rho, alpha, y, logOut)
   - eta1, eta2: REAL vectors length n
   - rho: REAL scalar or vector length n (recycled if length 1)
   - alpha: REAL matrix n×m (m>=2), column 1 is alpha1 (binary cut),
            columns 2..m are ordinal cuts (in increasing order per your R prep)
   - y: INT matrix n×2 (col1=y1∈{0,1}, col2=y2∈{0,1,2,...}; NOT normalized)
   - logOut: logical
*/
SEXP logLik_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha, SEXP Y, SEXP LogOut)
{
    /* coerce to expected types */
    Eta1  = PROTECT(coerceVector(Eta1,  REALSXP));
    Eta2  = PROTECT(coerceVector(Eta2,  REALSXP));
    Rho   = PROTECT(coerceVector(Rho,   REALSXP));
    Alpha = PROTECT(coerceVector(Alpha, REALSXP));
    Y     = PROTECT(coerceVector(Y,     INTSXP));
    LogOut= PROTECT(coerceVector(LogOut,LGLSXP));

    const double *eta1 = REAL(Eta1);
    const double *eta2 = REAL(Eta2);
    const double *rhoV = REAL(Rho);
    const double *A    = REAL(Alpha);
    const int    *Yp   = INTEGER(Y);

    /* dims */
    SEXP dimA = getAttrib(Alpha, R_DimSymbol);
    SEXP dimY = getAttrib(Y,     R_DimSymbol);
    if(TYPEOF(dimA) != INTSXP || TYPEOF(dimY) != INTSXP)
        error("alpha and y must have dim attributes");

    const int n   = INTEGER(dimA)[0];
    const int mA  = INTEGER(dimA)[1];   /* number of alpha columns (>=2) */
    const int nY  = INTEGER(dimY)[0];
    const int mY  = INTEGER(dimY)[1];

    if(nY != n || mY != 2) error("y must be an n x 2 matrix");
    if(mA < 2)             error("alpha must have at least two columns");
    if(XLENGTH(Eta1) != (R_xlen_t)n || XLENGTH(Eta2) != (R_xlen_t)n)
        error("eta1 and eta2 must have length n");

    const R_xlen_t nRho = XLENGTH(Rho);
    if(nRho != 1 && nRho != (R_xlen_t)n)
        warning("rho length (%ld) not equal to n (%d); recycling will be used",
                (long)nRho, n);

    const int logp = LOGICAL(LogOut)[0] == TRUE;
    const int steps = 128; /* Miwa grid; change here if you want another default */

    /* output */
    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *ll = REAL(out);

    for(int i = 0; i < n; ++i) {
        const int y1 = Yp[i + 0*n];   /* first column */
        const int y2 = Yp[i + 1*n];   /* second column (NOT normalized) */

        /* rho recycling exactly like R's vector indexing */
        const double rraw = rhoV[(nRho == 1) ? 0 : i];

        /* match mvtnorm: allow open interval (-1,1); clip slightly to avoid singularities */
        double r = rraw;
        if(r >=  1.0) r =  0.999999999;
        if(r <= -1.0) r = -0.999999999;

        /* A = alpha1 - eta1 */
        const double a1 = A[i + 0*n];
        const double Aup = a1 - eta1[i];

        if(y1 == 1) {
            /* P(Y1=1) = 1 - Phi(alpha1 - eta1) = upper tail */
            const double p = pnorm(Aup, 0.0, 1.0, /*lower.tail*/0, /*log.p*/0);
            ll[i] = logp ? ((p > 0.0) ? log(p) : R_NegInf) : (p > 0.0 ? p : 0.0);
            continue;
        }

        if(y1 != 0) {
            error("y1 must be 0 or 1 at row %d", i+1);
        }

        /* Build alpha2 vector: c(-Inf, alpha[i, 2:(mA)], +Inf)
           Indices in R: lower = alpha2[y2 + 1L], upper = alpha2[y2 + 2L]
           Here we compute B1 (lower) and B2 (upper) directly without materializing alpha2. */
        double B1, B2;

        /* lower bound */
        if(y2 < 0) {
            /* R would try to index alpha2[y2+1] with negative -> error; keep behavior */
            error("y2 < 0 at row %d", i+1);
        }
        if(y2 == 0) {
            B1 = R_NegInf;
        } else {
            /* y2 >= 1 -> alpha[i, y2+1] is the (y2+1)-th column (1-based) -> col index y2 (0-based) */
            const int colL = y2; /* 0-based: 0=alpha1, 1=alpha2, ... */
            if(colL >= mA)
                error("y2 too large for alpha at row %d (need column %d, have %d)", i+1, colL+1, mA);
            B1 = A[i + colL * n] - eta2[i];
        }

        /* upper bound */
        /* alpha2[y2+2]; if y2+2 == length(alpha2) -> +Inf */
        if(y2 + 1 >= mA) {
            /* top category: upper = +Inf */
            B2 = R_PosInf;
        } else {
            /* alpha[i, y2+2] -> 0-based column (y2+1) */
            const int colU = y2 + 1;
            if(colU >= mA)
                error("internal indexing error at row %d", i+1);
            B2 = A[i + colU * n] - eta2[i];
        }

        /* Probability: P(Z1 <= Aup, B1 < Z2 <= B2) = F(Aup, B2) - F(Aup, B1)
           with F(h,k) the BVN CDF up to (h,k) with corr r. */
        double p_up, p_lo;

        /* Short-cuts for infinite bounds (to match mvtnorm behavior) */
        if(!R_FINITE(B2) && B2 > 0) {
            /* Z2 <= +Inf -> F(Aup, +Inf) = Phi(Aup) */
            p_up = pnorm(Aup, 0.0, 1.0, 1, 0);
        } else {
            p_up = bvn_cdf_miwa(Aup, B2, r, steps);
        }

        if(!R_FINITE(B1) && B1 < 0) {
            /* Z2 <= -Inf -> 0 */
            p_lo = 0.0;
        } else {
            p_lo = bvn_cdf_miwa(Aup, B1, r, steps);
        }

        double p = p_up - p_lo;
        if(p < 0.0 && p > -1e-15) p = 0.0; /* guard tiny negatives */

        ll[i] = logp ? ((p > 0.0) ? log(p) : R_NegInf) : (p > 0.0 ? p : 0.0);
    }

    UNPROTECT(7);
    return out;
}

