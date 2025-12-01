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

/* --- helpers for rho link and BVN pdf ---------------------------------- */

/* d rho / d eta for rhogit inverse link, expressed in terms of rho */
static inline double rhogit_d1(double rho)
{
    double r = rho;
    if (r >=  1.0) r =  0.999999999;
    if (r <= -1.0) r = -0.999999999;
    double omr2 = 1.0 - r * r;
    if (omr2 <= 0.0) omr2 = DBL_MIN;
    /* (1 - rho^2)^(3/2) */
    return pow(omr2, 1.5);
}

/* bivariate standard normal pdf at (h, k) with correlation rho */
static inline double bvn_pdf(double h, double k, double rho)
{
    if (!R_FINITE(h) || !R_FINITE(k))
        return 0.0;

    double r = rho;
    if (r >=  1.0) r =  0.999999999;
    if (r <= -1.0) r = -0.999999999;
    double omr2 = 1.0 - r * r;
    if (omr2 <= 0.0) return 0.0;

    double quad = (h*h - 2.0*r*h*k + k*k) / omr2;
    double pref = 1.0 / (2.0 * M_PI * sqrt(omr2));

    return pref * exp(-0.5 * quad);
}

/* ---------------------------------------------------------------------- */
/* .Call interface:
 *
 *   score_rho_dontknow(Eta1, Eta2, Rho, Alpha, Y)
 *
 * Returns score wrt predictor eta_rho (length n vector).
 * - Eta1, Eta2: REAL vectors length n (means for latent Y1*, Y2*)
 * - Rho: REAL scalar or vector length n (parameter scale, after rhogit inv)
 * - Alpha: REAL matrix n x mA (mA>=2), col1=alpha1, cols2..mA=ordinal cuts
 * - Y: INT matrix n x 2 with columns (y1 in {0,1}, y2 in {0,1,2,...}) with y2
 *      already normalized to start at 0 (as in DK$pdf before calling C).
 */
SEXP score_rho_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha, SEXP Y)
{
    Eta1  = PROTECT(coerceVector(Eta1,  REALSXP));
    Eta2  = PROTECT(coerceVector(Eta2,  REALSXP));
    Rho   = PROTECT(coerceVector(Rho,   REALSXP));
    Alpha = PROTECT(coerceVector(Alpha, REALSXP));
    Y     = PROTECT(coerceVector(Y,     INTSXP));

    const double *eta1 = REAL(Eta1);
    const double *eta2 = REAL(Eta2);
    const double *rhoV = REAL(Rho);
    const double *A    = REAL(Alpha);
    const int    *Yp   = INTEGER(Y);

    SEXP dimA = getAttrib(Alpha, R_DimSymbol);
    SEXP dimY = getAttrib(Y,     R_DimSymbol);
    if (TYPEOF(dimA) != INTSXP || TYPEOF(dimY) != INTSXP)
        error("alpha and y must have dim attributes");

    const int n   = INTEGER(dimA)[0];
    const int mA  = INTEGER(dimA)[1];
    const int nY  = INTEGER(dimY)[0];
    const int mY  = INTEGER(dimY)[1];

    if (nY != n || mY != 2)
        error("y must be an n x 2 matrix");
    if (mA < 2)
        error("alpha must have at least two columns");
    if (XLENGTH(Eta1) != (R_xlen_t)n || XLENGTH(Eta2) != (R_xlen_t)n)
        error("eta1 and eta2 must have length n");

    const R_xlen_t nRho = XLENGTH(Rho);
    if (nRho != 1 && nRho != (R_xlen_t)n)
        warning("rho length (%ld) not equal to n (%d); recycling will be used",
                (long)nRho, n);

    const int steps = 128;

    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *sc = REAL(out);

    for (int i = 0; i < n; ++i) {
        const int y1 = Yp[i + 0*n];
        const int y2 = Yp[i + 1*n];

        if (y1 == 1) {
            /* dk: does not depend on rho */
            sc[i] = 0.0;
            continue;
        }
        if (y1 != 0) {
            error("y1 must be 0 or 1 at row %d", i+1);
        }

        /* rho for this obs (with recycling) */
        const double rraw = rhoV[(nRho == 1) ? 0 : i];
        double r = rraw;
        if (r >=  1.0) r =  0.999999999;
        if (r <= -1.0) r = -0.999999999;

        /* A = alpha1 - eta1 */
        const double a1 = A[i + 0*n];
        const double Aup = a1 - eta1[i];

        /* build B1, B2 from observed y2 (same logic as in logLik_dontknow) */
        double B1, B2;

        if (y2 < 0) {
            error("y2 < 0 at row %d", i+1);
        }
        if (y2 == 0) {
            B1 = R_NegInf;
        } else {
            const int colL = y2;   /* 0-based: 0=alpha1, 1=alpha2, ... */
            if (colL >= mA)
                error("y2 too large for alpha at row %d (need column %d, have %d)",
                      i+1, colL+1, mA);
            B1 = A[i + colL * n] - eta2[i];
        }

        if (y2 + 1 >= mA) {
            B2 = R_PosInf;
        } else {
            const int colU = y2 + 1;
            if (colU >= mA)
                error("internal indexing error at row %d", i+1);
            B2 = A[i + colU * n] - eta2[i];
        }

        /* rectangle probability: p = F(Aup,B2) - F(Aup,B1) */
        double p_up, p_lo;
        if (!R_FINITE(B2) && B2 > 0.0) {
            p_up = pnorm(Aup, 0.0, 1.0, 1, 0);
        } else {
            p_up = bvn_cdf_miwa(Aup, B2, r, steps);
        }
        if (!R_FINITE(B1) && B1 < 0.0) {
            p_lo = 0.0;
        } else {
            p_lo = bvn_cdf_miwa(Aup, B1, r, steps);
        }

        double p = p_up - p_lo;
        if (p < 0.0 && p > -1e-15)
            p = 0.0;

        if (p <= 0.0 || !R_FINITE(p)) {
            sc[i] = 0.0;
        } else {
            /* derivative dp/d rho = phi2(Aup,B2;r) - phi2(Aup,B1;r) */
            double phi_up = bvn_pdf(Aup, B2, r);
            double phi_lo = bvn_pdf(Aup, B1, r);
            double p_prime = phi_up - phi_lo;

            double s_rho = p_prime / p;  /* d l / d rho */
            double d1    = rhogit_d1(r); /* d rho / d eta_rho */

            sc[i] = s_rho * d1;          /* d l / d eta_rho */
        }
    }

    UNPROTECT(6);
    return out;
}

/* ---------------------------------------------------------------------- */
/* .Call interface:
 *
 *   hess_rho_dontknow(Eta1, Eta2, Rho, Alpha)
 *
 * Returns Fisher information (expected negative Hessian) wrt predictor
 * eta_rho, i.e. per observation weight > 0.
 *
 * For each obs i:
 *   I_rho[i]  = sum_c (p'_c^2 / p_c),  over c=0..K-1
 *   I_eta[i]  = (d rho / d eta)^2 * I_rho[i]
 */
SEXP hess_rho_dontknow(SEXP Eta1, SEXP Eta2, SEXP Rho, SEXP Alpha)
{
    Eta1  = PROTECT(coerceVector(Eta1,  REALSXP));
    Eta2  = PROTECT(coerceVector(Eta2,  REALSXP));
    Rho   = PROTECT(coerceVector(Rho,   REALSXP));
    Alpha = PROTECT(coerceVector(Alpha, REALSXP));

    const double *eta1 = REAL(Eta1);
    const double *eta2 = REAL(Eta2);
    const double *rhoV = REAL(Rho);
    const double *A    = REAL(Alpha);

    SEXP dimA = getAttrib(Alpha, R_DimSymbol);
    if (TYPEOF(dimA) != INTSXP)
        error("alpha must have dim attributes");

    const int n   = INTEGER(dimA)[0];
    const int mA  = INTEGER(dimA)[1];
    if (mA < 2)
        error("alpha must have at least two columns");
    if (XLENGTH(Eta1) != (R_xlen_t)n || XLENGTH(Eta2) != (R_xlen_t)n)
        error("eta1 and eta2 must have length n");

    const R_xlen_t nRho = XLENGTH(Rho);
    if (nRho != 1 && nRho != (R_xlen_t)n)
        warning("rho length (%ld) not equal to n (%d); recycling will be used",
                (long)nRho, n);

    const int K     = mA - 1;  /* number of ordinal categories */
    const int steps = 128;

    SEXP out = PROTECT(allocVector(REALSXP, n));
    double *w = REAL(out);

    for (int i = 0; i < n; ++i) {
        const double rraw = rhoV[(nRho == 1) ? 0 : i];
        double r = rraw;
        if (r >=  1.0) r =  0.999999999;
        if (r <= -1.0) r = -0.999999999;

        /* A = alpha1 - eta1 */
        const double a1 = A[i + 0*n];
        const double Aup = a1 - eta1[i];

        double I_rho = 0.0;

        /* loop over ordinal categories c = 0..K-1 */
        for (int c = 0; c < K; ++c) {
            double B1, B2;

            /* replicate alpha2 = c(-Inf, alpha2..alphaK, +Inf) indexing
             * B1 = alpha2[c+1] - eta2, B2 = alpha2[c+2] - eta2
             */
            if (c == 0) {
                B1 = R_NegInf;
            } else {
                /* alpha2[c+1] = alpha[i, c+1] (1-based), 0-based col index c */
                if (c >= mA)
                    error("internal indexing error in hess_rho_dontknow (B1) at row %d", i+1);
                B1 = A[i + c * n] - eta2[i];
            }

            if (c + 1 >= K) {
                /* top category -> +Inf */
                B2 = R_PosInf;
            } else {
                /* alpha2[c+2] = alpha[i, c+2] -> 0-based col index (c+1) */
                int colU = c + 1;
                if (colU >= mA)
                    error("internal indexing error in hess_rho_dontknow (B2) at row %d", i+1);
                B2 = A[i + colU * n] - eta2[i];
            }

            /* p_c = F(Aup,B2) - F(Aup,B1) */
            double p_up, p_lo;
            if (!R_FINITE(B2) && B2 > 0.0) {
                p_up = pnorm(Aup, 0.0, 1.0, 1, 0);
            } else {
                p_up = bvn_cdf_miwa(Aup, B2, r, steps);
            }
            if (!R_FINITE(B1) && B1 < 0.0) {
                p_lo = 0.0;
            } else {
                p_lo = bvn_cdf_miwa(Aup, B1, r, steps);
            }

            double p_c = p_up - p_lo;
            if (p_c < 0.0 && p_c > -1e-15)
                p_c = 0.0;

            if (p_c <= 0.0 || !R_FINITE(p_c))
                continue;

            /* p'_c = phi2(Aup,B2;r) - phi2(Aup,B1;r) */
            double phi_up = bvn_pdf(Aup, B2, r);
            double phi_lo = bvn_pdf(Aup, B1, r);
            double p_prime_c = phi_up - phi_lo;

            I_rho += (p_prime_c * p_prime_c) / p_c;
        }

        /* transform to eta_rho-scale: I_eta = (d rho / d eta)^2 * I_rho */
        double d1 = rhogit_d1(r);
        double w_i = d1 * d1 * I_rho;
        if (!R_FINITE(w_i) || w_i <= 0.0)
            w_i = 1e-6;

        w[i] = w_i;
    }

    UNPROTECT(5);
    return out;
}
