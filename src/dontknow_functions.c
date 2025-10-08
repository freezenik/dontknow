#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
 * ===================================================================================
 * Bivariate Normal CDF
 *
 * C translation of Alan Genz's `bvn.f` Fortran routine. This is used by the
 * R package `mvtnorm` for bivariate normal probabilities, so using this
 * implementation ensures the results will be virtually identical.
 *
 * Original Fortran code available from Alan Genz's website:
 * http://www.math.wsu.edu/faculty/genz/software/software.html
 * ===================================================================================
 */
static double bvn(double h, double k, double r) {
    if (r == 0) {
        return pnorm(h, 0.0, 1.0, 1, 0) * pnorm(k, 0.0, 1.0, 1, 0);
    }

    double hk = h * k;
    double bvn = 0.0;

    if (fabs(r) < 0.925) {
        double hks = (h * h + k * k) / 2.0;
        double r2 = r * r;
        if (fabs(r) > 0.75) { // series for |r| > 0.75
            double g = (1.0 - r2);
            bvn = exp(-(hks - hk * r) / g) / (2.0 * M_PI * sqrt(g));
            double hr = (h - k * r) / sqrt(g);
            double kr = (k - h * r) / sqrt(g);
            if (hr > 0) bvn *= pnorm(hr, 0.0, 1.0, 0, 0);
            if (kr > 0) bvn *= pnorm(kr, 0.0, 1.0, 0, 0);
        } else { // series for |r| < 0.75
            bvn = exp(-hks) * (1.0 - r * r * (1.0 - hk * r / 3.0) +
                        r2 * r2 * (0.5 - 0.25 * r2 + hk * r * (1.0 - 0.25 * r2) / 3.0));
        }
    } else { // direct integration for |r| > 0.925
        if (r < 0) {
            k = -k;
            hk = -hk;
        }
        if (fabs(h - k) < 0.001 * (fabs(h) + fabs(k))) {
            double s = (h + k) / 2.0;
            double d = (h - k) / 2.0;
            double hs = s / sqrt(2.0 * (1.0 + r));
            bvn = pnorm(hs, 0.0, 1.0, 1, 0) -
                  exp(-(s * s) / (1.0 + r)) * (d * d / (1.0 - r) + 1.0) * 0.05641895835;
        } else {
             bvn = pnorm(h, 0.0, 1.0, 1, 0) * pnorm(k, 0.0, 1.0, 1, 0) -
                  exp(-(hk) / (2.0 * r)) * hk * (1.0 - hk / (4.0 * r)) / (24.0 * r);
        }
    }

    return bvn;
}


/*
 * ===================================================================================
 * Main log-likelihood function to be called from R.
 *
 * This function calculates the log-likelihood for a bivariate model where one
 * outcome is binary (Y1) and the other is ordinal (Y2). The case where Y1=1
 * corresponds to "Don't Know".
 *
 * Corresponds to the R function `logLik_dontknow()`.
 * ===================================================================================
 */
SEXP logLik_dontknow(SEXP s_eta1, SEXP s_eta2, SEXP s_rho, SEXP s_alpha, SEXP s_y, SEXP s_log)
{
    // --- Argument Coercion and Pointer Setup ---
    s_eta1 = PROTECT(coerceVector(s_eta1, REALSXP));
    s_eta2 = PROTECT(coerceVector(s_eta2, REALSXP));
    s_rho  = PROTECT(coerceVector(s_rho, REALSXP));
    s_alpha = PROTECT(coerceVector(s_alpha, REALSXP));
    s_y = PROTECT(coerceVector(s_y, INTSXP)); // y is integer
    s_log = PROTECT(coerceVector(s_log, LGLSXP));

    double *eta1  = REAL(s_eta1);
    double *eta2  = REAL(s_eta2);
    double *rho   = REAL(s_rho);
    double *alpha = REAL(s_alpha);
    int *y = INTEGER(s_y);
    int log_p = LOGICAL(s_log)[0];

    // --- Get Dimensions ---
    int n = length(s_eta1);
    int n_alpha_cols = INTEGER(getAttrib(s_alpha, R_DimSymbol))[1];

    // --- Allocate Result Vector ---
    SEXP s_ll = PROTECT(allocVector(REALSXP, n));
    double *ll = REAL(s_ll);

    // --- Main Loop ---
    for (int i = 0; i < n; ++i) {
        int y1 = y[i];

        if (y1 == 1) {
            // Case 1: Y1 = 1 ("Don't Know")
            // P(Y1 = 1) = 1 - Phi(alpha1 - eta1) = P(Z1 > alpha1 - eta1)
            // In R: pnorm(alpha[i, 1] - eta1[i], lower.tail = FALSE)
            // Matrix access: alpha is column-major. alpha[i, 1] is alpha[i + 0*n].
            double val = alpha[i] - eta1[i];
            ll[i] = pnorm(val, 0.0, 1.0, 0, log_p); // lower.tail=0, log.p=log_p
        } else {
            // Case 2: Y1 = 0
            // Here we compute P(Z1 < A, B1 < Z2 < B2)
            // = P(Z1 < A, Z2 < B2) - P(Z1 < A, Z2 < B1)
            int y2 = y[i + n]; // y is a matrix, y[,2] starts at index n
            double current_rho = rho[i];

            // Upper bound for the first dimension (Z1)
            double A = alpha[i] - eta1[i];

            // Bounds for the second dimension (Z2)
            double B1, B2;

            // Lower bound B1 from alpha_{2, y2}
            if (y2 == 0) {
                B1 = R_NegInf;
            } else {
                // R code uses alpha2[y2[i] + 1L].
                // alpha2 is c(-Inf, alpha[i, 2], alpha[i, 3], ...).
                // So for y2=1, index is 2 -> alpha[i, 2].
                // alpha column index is y2 (0-indexed) + 1 (for R), -1 (for C) -> y2.
                B1 = alpha[i + y2 * n] - eta2[i];
            }

            // Upper bound B2 from alpha_{2, y2+1}
            if (y2 == n_alpha_cols - 1) { // y2 is 0-indexed, highest category.
                B2 = R_PosInf;
            } else {
                // R code uses alpha2[y2[i] + 2L].
                // For y2=0, index is 2 -> alpha[i, 2] -> C col index 1.
                B2 = alpha[i + (y2 + 1) * n] - eta2[i];
            }

            // --- Calculate Probabilities ---
            double p_up, p_lo;

            // P(Z1 < A, Z2 < B2)
            if (B2 == R_PosInf) {
                p_up = pnorm(A, 0.0, 1.0, 1, 0); // P(Z1 < A)
            } else {
                p_up = bvn(A, B2, current_rho);
            }

            // P(Z1 < A, Z2 < B1)
            if (B1 == R_NegInf) {
                p_lo = 0.0;
            } else {
                p_lo = bvn(A, B1, current_rho);
            }

            double p = p_up - p_lo;

            // --- Store result with guards ---
            if (log_p) {
                ll[i] = (p > 0) ? log(p) : R_NegInf;
            } else {
                ll[i] = (p > 0) ? p : 0.0;
            }
        }
    }

    UNPROTECT(7);
    return s_ll;
}

